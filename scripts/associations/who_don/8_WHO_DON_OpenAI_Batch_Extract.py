#!/usr/bin/env python3
"""
Prepare and run OpenAI Batch workflows for WHO DON country extraction or
adjudication, then parse the resulting country candidates back into repo
artifacts.

Default workflow:
  1. prepare  -> build Batch API JSONL requests from WHO DON input
  2. submit   -> upload the request file and create a batch job
  3. status   -> inspect batch progress
  4. download -> fetch output/error files for a completed batch
  5. parse    -> validate and flatten model outputs into CSV/JSONL artifacts

Required environment variable for API calls:
  OPENAI_API_KEY
"""

from __future__ import annotations

import argparse
import csv
import json
import mimetypes
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "pathogen_association_data" / "WHO" / "disease_outbreak_news"
BATCH_DIR = OUTPUT_DIR / "openai_batch"
INPUT_JSONL = OUTPUT_DIR / "who_don_unresolved_llm_input.jsonl"
COUNTRY_ALIAS_CSV = Path(__file__).resolve().parent / "who_don_country_aliases.csv"

DEFAULT_MODEL = "gpt-5.4-mini-2026-03-17"
DEFAULT_COMPLETION_WINDOW = "24h"
DEFAULT_ENDPOINT = "/v1/responses"
API_BASE = "https://api.openai.com/v1"


def ensure_batch_dir() -> Path:
    BATCH_DIR.mkdir(parents=True, exist_ok=True)
    return BATCH_DIR


def message(*parts: Any) -> None:
    print(*parts, flush=True)


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSONL at {path}:{line_number}: {exc}") from exc
    return rows


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False))
            handle.write("\n")


def load_dotenv_values(dotenv_path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not dotenv_path.exists():
        return values

    for raw_line in dotenv_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue

        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

        if not key:
            continue

        if "#" in value and not (value.startswith('"') or value.startswith("'")):
            value = value.split("#", 1)[0].strip()

        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]

        values[key] = value

    return values


def load_country_allowlist() -> list[str]:
    countries: set[str] = set()
    with COUNTRY_ALIAS_CSV.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            geography_type = (row.get("geography_type") or "").strip()
            country_standard = (row.get("country_standard") or "").strip()
            is_ambiguous = (row.get("is_ambiguous") or "").strip().lower()
            if geography_type == "country" and country_standard and is_ambiguous != "true":
                countries.add(country_standard)

    base_countries = [
        "Afghanistan", "Albania", "Algeria", "Andorra", "Angola", "Antigua and Barbuda",
        "Argentina", "Armenia", "Australia", "Austria", "Azerbaijan", "Bahamas", "Bahrain",
        "Bangladesh", "Barbados", "Belarus", "Belgium", "Belize", "Benin", "Bhutan", "Bolivia",
        "Bosnia and Herzegovina", "Botswana", "Brazil", "Brunei", "Bulgaria", "Burkina Faso",
        "Burundi", "Cabo Verde", "Cambodia", "Cameroon", "Canada", "Central African Republic",
        "Chad", "Chile", "China", "Colombia", "Comoros", "Costa Rica", "Croatia", "Cuba",
        "Cyprus", "Czechia", "Denmark", "Djibouti", "Dominica", "Dominican Republic", "Ecuador",
        "Egypt", "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia", "Eswatini",
        "Ethiopia", "Fiji", "Finland", "France", "Gabon", "Gambia", "Georgia", "Germany",
        "Ghana", "Greece", "Grenada", "Guatemala", "Guinea", "Guinea-Bissau", "Guyana",
        "Haiti", "Honduras", "Hungary", "Iceland", "India", "Indonesia", "Iran", "Iraq",
        "Ireland", "Israel", "Italy", "Jamaica", "Japan", "Jordan", "Kazakhstan", "Kenya",
        "Kiribati", "Kuwait", "Kyrgyzstan", "Laos", "Latvia", "Lebanon", "Lesotho", "Liberia",
        "Libya", "Liechtenstein", "Lithuania", "Luxembourg", "Madagascar", "Malawi", "Malaysia",
        "Maldives", "Mali", "Malta", "Marshall Islands", "Mauritania", "Mauritius", "Mexico",
        "Micronesia", "Moldova", "Monaco", "Mongolia", "Montenegro", "Morocco", "Mozambique",
        "Myanmar", "Namibia", "Nauru", "Nepal", "Netherlands", "New Zealand", "Nicaragua",
        "Niger", "Nigeria", "North Korea", "North Macedonia", "Norway", "Oman", "Pakistan",
        "Palau", "Panama", "Papua New Guinea", "Paraguay", "Peru", "Philippines", "Poland",
        "Portugal", "Qatar", "Republic of the Congo", "Romania", "Russia", "Rwanda",
        "Saint Kitts and Nevis", "Saint Lucia", "Saint Vincent and the Grenadines", "Samoa",
        "San Marino", "Sao Tome and Principe", "Saudi Arabia", "Senegal", "Serbia",
        "Seychelles", "Sierra Leone", "Singapore", "Slovakia", "Slovenia", "Solomon Islands",
        "Somalia", "South Africa", "South Korea", "South Sudan", "Spain", "Sri Lanka", "Sudan",
        "Suriname", "Sweden", "Switzerland", "Syria", "Tajikistan", "Tanzania", "Thailand",
        "Timor-Leste", "Togo", "Tonga", "Trinidad and Tobago", "Tunisia", "Turkey",
        "Turkmenistan", "Tuvalu", "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom",
        "United States", "Uruguay", "Uzbekistan", "Vanuatu", "Venezuela", "Vietnam", "Yemen",
        "Zambia", "Zimbabwe", "Palestine", "Holy See", "Kosovo", "Democratic Republic of the Congo",
        "Cote d'Ivoire"
    ]
    countries.update(base_countries)
    return sorted(countries)


def extraction_schema(country_allowlist: list[str]) -> dict[str, Any]:
    return {
        "name": "who_don_country_extraction",
        "strict": True,
        "schema": {
            "type": "object",
            "additionalProperties": False,
            "required": [
                "record_key",
                "has_country_evidence",
                "country_evidence",
                "reasoning_label",
                "confidence",
            ],
            "properties": {
                "record_key": {"type": "string"},
                "has_country_evidence": {"type": "boolean"},
                "country_evidence": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "additionalProperties": False,
                        "required": ["country", "evidence_span"],
                        "properties": {
                            "country": {"type": "string", "enum": country_allowlist},
                            "evidence_span": {"type": "string"},
                        },
                    },
                },
                "reasoning_label": {
                    "type": "string",
                    "enum": [
                        "explicit_event_country",
                        "background_only",
                        "regional_only",
                        "no_country",
                    ],
                },
                "confidence": {
                    "type": "string",
                    "enum": ["high", "medium", "low"],
                },
            },
        },
    }


def system_prompt(workflow: str) -> str:
    if workflow == "adjudicate":
        return (
            "You validate country evidence from WHO Disease Outbreak News records. "
            "Return only structured JSON matching the schema. "
            "This is a strict adjudication pass, not a recall-maximizing extraction pass. "
            "Keep only countries explicitly supported by the record text itself. "
            "Candidate countries or spans shown in the prompt may be wrong; use them only as hints. "
            "Return at most one strongest evidence_span per retained country. "
            "Prefer short verbatim spans that explicitly contain the country name. "
            "If an explicit country-literal span exists, do not use weaker travel-only, background-only, "
            "or indirect location spans instead. "
            "Do not infer countries from disease history, travel assumptions, subnational geography knowledge, "
            "or regional labels. Do not turn regions such as West Africa or Region of the Americas into country lists. "
            "If the support is weak or only background context, drop the country or return no_country. "
            "If no explicit country evidence is present, return has_country_evidence=false, "
            "country_evidence=[], reasoning_label='no_country'. "
            "Evidence spans must be short verbatim snippets from the record."
        )

    return (
        "You extract country evidence from WHO Disease Outbreak News records. "
        "Return only structured JSON matching the schema. "
        "Use only countries explicitly stated in the record text. "
        "Do not infer countries from disease history, travel assumptions, or regional labels. "
        "Do not turn regions such as West Africa or Region of the Americas into country lists. "
        "Return one object per country in country_evidence, with keys country and evidence_span. "
        "If no explicit country evidence is present, return has_country_evidence=false, "
        "country_evidence=[], reasoning_label='no_country'. "
        "Evidence spans must be short verbatim snippets from the record."
    )


def render_prompt_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, indent=2)
    return str(value)


def build_user_prompt(payload: dict[str, Any], workflow: str) -> str:
    input_obj = payload.get("input", {})
    parts = [
        f"record_key: {payload.get('record_key', '')}",
        f"Title: {input_obj.get('Title', '')}",
        f"article_url: {input_obj.get('article_url', '')}",
        "",
        "summary_text:",
        input_obj.get("summary_text", "") or "",
        "",
        "overview_text:",
        input_obj.get("overview_text", "") or "",
        "",
        "response_text:",
        input_obj.get("response_text", "") or "",
        "",
        "slug_hint:",
        input_obj.get("slug_hint", "") or "",
        "",
        "region_hint:",
        input_obj.get("region_hint", "") or "",
    ]

    optional_fields = [
        "current_best_country_list",
        "current_best_country_rows",
        "manual_review_country_rows",
        "adjudication_reason",
        "adjudication_notes",
        "llm_country_count",
    ]

    if workflow == "adjudicate":
        for field_name in optional_fields:
            if field_name in input_obj:
                parts.extend(
                    [
                        "",
                        f"{field_name}:",
                        render_prompt_value(input_obj.get(field_name)),
                    ]
                )

    return "\n".join(parts).strip()


def build_responses_body(
    payload: dict[str, Any],
    model: str,
    country_allowlist: list[str],
    reasoning_effort: str | None,
    workflow: str,
) -> dict[str, Any]:
    body: dict[str, Any] = {
        "model": model,
        "instructions": system_prompt(workflow),
        "input": build_user_prompt(payload, workflow=workflow),
        "text": {
            "format": {
                "type": "json_schema",
                **extraction_schema(country_allowlist),
            }
        },
    }
    if reasoning_effort and reasoning_effort != "none":
        body["reasoning"] = {"effort": reasoning_effort}
    return body


def build_chat_completions_body(
    payload: dict[str, Any],
    model: str,
    country_allowlist: list[str],
    reasoning_effort: str | None,
    workflow: str,
) -> dict[str, Any]:
    body: dict[str, Any] = {
        "model": model,
        "temperature": 0,
        "messages": [
            {"role": "system", "content": system_prompt(workflow)},
            {"role": "user", "content": build_user_prompt(payload, workflow=workflow)},
        ],
        "response_format": {
            "type": "json_schema",
            "json_schema": extraction_schema(country_allowlist),
        },
    }
    if reasoning_effort and reasoning_effort != "none":
        body["reasoning"] = {"effort": reasoning_effort}
    return body


def batch_request_row(
    payload: dict[str, Any],
    model: str,
    country_allowlist: list[str],
    reasoning_effort: str | None,
    endpoint: str,
    workflow: str,
) -> dict[str, Any]:
    if endpoint == "/v1/responses":
        body = build_responses_body(
            payload=payload,
            model=model,
            country_allowlist=country_allowlist,
            reasoning_effort=reasoning_effort,
            workflow=workflow,
        )
    elif endpoint == "/v1/chat/completions":
        body = build_chat_completions_body(
            payload=payload,
            model=model,
            country_allowlist=country_allowlist,
            reasoning_effort=reasoning_effort,
            workflow=workflow,
        )
    else:
        raise ValueError(f"Unsupported endpoint for batch requests: {endpoint}")

    return {
        "custom_id": str(payload.get("record_key", "")),
        "method": "POST",
        "url": endpoint,
        "body": body,
    }


def load_api_key() -> str:
    env_candidates = [
        os.getenv("OPENAI_API_KEY", "").strip(),
        os.getenv("openai_api_key", "").strip(),
    ]
    api_key = next((value for value in env_candidates if value), "")

    if not api_key:
        dotenv_values = load_dotenv_values(REPO_ROOT / ".env")
        dotenv_candidates = [
            dotenv_values.get("OPENAI_API_KEY", "").strip(),
            dotenv_values.get("openai_api_key", "").strip(),
        ]
        api_key = next((value for value in dotenv_candidates if value), "")

    if not api_key:
        raise RuntimeError(
            "OpenAI API key is required for submit/status/download actions. "
            "Set OPENAI_API_KEY or openai_api_key in the environment or .env."
        )
    return api_key


def api_request(
    method: str,
    path: str,
    api_key: str,
    json_body: dict[str, Any] | None = None,
    content_type: str = "application/json",
    body_bytes: bytes | None = None,
) -> dict[str, Any]:
    url = f"{API_BASE}{path}"
    if json_body is not None:
        body_bytes = json.dumps(json_body).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
    }
    if body_bytes is not None:
        headers["Content-Type"] = content_type
    request = urllib.request.Request(url, data=body_bytes, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as response:
            raw = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        raw = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"OpenAI API error {exc.code} for {path}: {raw}") from exc
    return json.loads(raw)


def multipart_form_data(fields: dict[str, str], files: list[tuple[str, Path, str]]) -> tuple[bytes, str]:
    boundary = f"----CodexBoundary{uuid.uuid4().hex}"
    chunks: list[bytes] = []

    for name, value in fields.items():
        chunks.extend(
            [
                f"--{boundary}\r\n".encode("utf-8"),
                f'Content-Disposition: form-data; name="{name}"\r\n\r\n'.encode("utf-8"),
                value.encode("utf-8"),
                b"\r\n",
            ]
        )

    for field_name, file_path, mime_type in files:
        filename = file_path.name
        file_bytes = file_path.read_bytes()
        chunks.extend(
            [
                f"--{boundary}\r\n".encode("utf-8"),
                (
                    f'Content-Disposition: form-data; name="{field_name}"; '
                    f'filename="{filename}"\r\n'
                ).encode("utf-8"),
                f"Content-Type: {mime_type}\r\n\r\n".encode("utf-8"),
                file_bytes,
                b"\r\n",
            ]
        )

    chunks.append(f"--{boundary}--\r\n".encode("utf-8"))
    return b"".join(chunks), f"multipart/form-data; boundary={boundary}"


def upload_file(file_path: Path, api_key: str, purpose: str = "batch") -> dict[str, Any]:
    mime_type = mimetypes.guess_type(file_path.name)[0] or "application/octet-stream"
    body, content_type = multipart_form_data(
        fields={"purpose": purpose},
        files=[("file", file_path, mime_type)],
    )
    return api_request(
        method="POST",
        path="/files",
        api_key=api_key,
        body_bytes=body,
        content_type=content_type,
    )


def download_file_content(file_id: str, api_key: str) -> bytes:
    url = f"{API_BASE}/files/{file_id}/content"
    request = urllib.request.Request(
        url,
        headers={"Authorization": f"Bearer {api_key}"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(request) as response:
            return response.read()
    except urllib.error.HTTPError as exc:
        raw = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"OpenAI API error {exc.code} for file {file_id}: {raw}") from exc


def prepare_requests(args: argparse.Namespace) -> None:
    ensure_batch_dir()
    input_path = Path(args.input_jsonl)
    rows = read_jsonl(input_path)

    if args.record_keys:
        wanted = {key.strip() for key in args.record_keys.split(",") if key.strip()}
        rows = [row for row in rows if str(row.get("record_key", "")).strip() in wanted]

    if args.limit is not None:
        rows = rows[: args.limit]

    if not rows:
        raise ValueError("No input rows remain after applying --record-keys/--limit filters.")

    country_allowlist = load_country_allowlist()

    batch_rows = [
        batch_request_row(
            payload=row,
            model=args.model,
            country_allowlist=country_allowlist,
            reasoning_effort=args.reasoning_effort,
            endpoint=args.endpoint,
            workflow=args.workflow,
        )
        for row in rows
    ]

    request_path = Path(args.requests_jsonl)
    write_jsonl(request_path, batch_rows)

    manifest = {
        "prepared_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_jsonl": str(input_path),
        "requests_jsonl": str(request_path),
        "model": args.model,
        "endpoint": args.endpoint,
        "workflow": args.workflow,
        "reasoning_effort": args.reasoning_effort,
        "n_requests": len(batch_rows),
        "limit": args.limit,
        "record_keys": args.record_keys,
    }
    manifest_path = request_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    message(f"Prepared OpenAI batch requests: {request_path}")
    message(f"Requests written: {len(batch_rows)}")


def submit_batch(args: argparse.Namespace) -> None:
    api_key = load_api_key()
    ensure_batch_dir()
    requests_path = Path(args.requests_jsonl)
    upload_result = upload_file(requests_path, api_key=api_key, purpose="batch")
    batch_result = api_request(
        method="POST",
        path="/batches",
        api_key=api_key,
        json_body={
            "input_file_id": upload_result["id"],
            "endpoint": args.endpoint,
            "completion_window": args.completion_window,
            "metadata": {
                "workflow": (
                    "who_don_country_adjudication"
                    if args.workflow == "adjudicate"
                    else "who_don_country_extraction"
                ),
                "requests_file": requests_path.name,
                "model": args.model,
                "endpoint": args.endpoint,
                "prompt_workflow": args.workflow,
            },
        },
    )

    upload_path = Path(args.upload_json)
    batch_path = Path(args.batch_json)
    upload_path.write_text(json.dumps(upload_result, indent=2), encoding="utf-8")
    batch_path.write_text(json.dumps(batch_result, indent=2), encoding="utf-8")

    message(f"Uploaded batch input file: {upload_result['id']}")
    message(f"Created batch job: {batch_result['id']}")
    message(f"Batch metadata saved to: {batch_path}")


def fetch_batch_status(args: argparse.Namespace) -> None:
    api_key = load_api_key()
    batch = api_request("GET", f"/batches/{args.batch_id}", api_key=api_key)
    output_path = Path(args.batch_json)
    output_path.write_text(json.dumps(batch, indent=2), encoding="utf-8")
    message(json.dumps(batch, indent=2))


def download_batch_outputs(args: argparse.Namespace) -> None:
    api_key = load_api_key()
    batch_info = json.loads(Path(args.batch_json).read_text(encoding="utf-8"))
    output_file_id = batch_info.get("output_file_id")
    error_file_id = batch_info.get("error_file_id")

    if output_file_id:
        output_bytes = download_file_content(output_file_id, api_key=api_key)
        Path(args.output_jsonl).write_bytes(output_bytes)
        message(f"Saved batch output file to: {args.output_jsonl}")
    else:
        message("No output_file_id present yet.")

    if error_file_id:
        error_bytes = download_file_content(error_file_id, api_key=api_key)
        Path(args.error_jsonl).write_bytes(error_bytes)
        message(f"Saved batch error file to: {args.error_jsonl}")
    else:
        message("No error_file_id present.")


def extract_responses_output_content(response_body: dict[str, Any]) -> str:
    outputs = response_body.get("output", [])
    if not isinstance(outputs, list) or not outputs:
        raise ValueError("No output array in Responses API batch response.")

    text_parts: list[str] = []
    refusals: list[str] = []

    for output in outputs:
        if output.get("type") != "message":
            continue
        for item in output.get("content", []):
            item_type = item.get("type")
            if item_type == "output_text" and "text" in item:
                text_parts.append(item["text"])
            elif item_type == "refusal" and "refusal" in item:
                refusals.append(item["refusal"])

    if text_parts:
        return "".join(text_parts)
    if refusals:
        raise ValueError(f"Model refusal: {' | '.join(refusals)}")
    raise ValueError("Could not extract output_text from Responses API response body.")


def extract_chat_completions_content(response_body: dict[str, Any]) -> str:
    choices = response_body.get("choices", [])
    if not choices:
        raise ValueError("No choices in batch response body.")

    message_obj = choices[0].get("message", {})
    content = message_obj.get("content")

    if isinstance(content, str):
        return content

    if isinstance(content, list):
        text_parts = []
        for item in content:
            if isinstance(item, dict):
                if item.get("type") in {"output_text", "text"} and "text" in item:
                    text_parts.append(item["text"])
        if text_parts:
            return "".join(text_parts)

    if message_obj.get("refusal"):
        raise ValueError(f"Model refusal: {message_obj['refusal']}")

    raise ValueError("Could not extract text content from Chat Completions response body.")


def extract_model_output_content(response_body: dict[str, Any]) -> str:
    if "output" in response_body:
        return extract_responses_output_content(response_body)
    if "choices" in response_body:
        return extract_chat_completions_content(response_body)
    raise ValueError("Unknown batch response body shape.")


def validate_extraction_row(
    row: dict[str, Any],
    country_allowlist: set[str],
) -> tuple[bool, str]:
    required = {
        "record_key",
        "has_country_evidence",
        "country_evidence",
        "reasoning_label",
        "confidence",
    }
    missing = required.difference(row.keys())
    if missing:
        return False, f"Missing required keys: {sorted(missing)}"

    country_evidence = row.get("country_evidence", [])
    if not isinstance(country_evidence, list):
        return False, "country_evidence must be an array."
    for item in country_evidence:
        if not isinstance(item, dict):
            return False, "country_evidence entries must be objects."
        if set(item.keys()) != {"country", "evidence_span"}:
            return False, "country_evidence entries must contain only country and evidence_span."
        if item["country"] not in country_allowlist:
            return False, "One or more countries are outside the allowlist."
        if not isinstance(item["evidence_span"], str):
            return False, "evidence_span must be a string."
    if row["has_country_evidence"] is False and country_evidence:
        return False, "has_country_evidence=false but country_evidence is present."
    if row["reasoning_label"] == "no_country" and country_evidence:
        return False, "reasoning_label=no_country but country_evidence is present."
    return True, ""


def parse_batch_outputs(args: argparse.Namespace) -> None:
    ensure_batch_dir()
    input_rows = {row.get("record_key"): row for row in read_jsonl(Path(args.input_jsonl))}
    output_rows = read_jsonl(Path(args.output_jsonl))
    country_allowlist = set(load_country_allowlist())

    parsed_records: list[dict[str, Any]] = []
    flattened_rows: list[dict[str, Any]] = []
    invalid_rows: list[dict[str, Any]] = []

    for batch_row in output_rows:
        custom_id = batch_row.get("custom_id")
        response = batch_row.get("response") or {}
        error = batch_row.get("error")

        if error:
            invalid_rows.append(
                {
                    "record_key": custom_id,
                    "error_type": "batch_error",
                    "error_detail": json.dumps(error, ensure_ascii=False),
                }
            )
            continue

        body = response.get("body", {})
        try:
            content = extract_model_output_content(body)
            parsed = json.loads(content)
        except Exception as exc:  # noqa: BLE001
            invalid_rows.append(
                {
                    "record_key": custom_id,
                    "error_type": "parse_error",
                    "error_detail": str(exc),
                }
            )
            continue

        is_valid, validation_message = validate_extraction_row(parsed, country_allowlist)
        if not is_valid:
            invalid_rows.append(
                {
                    "record_key": custom_id,
                    "error_type": "validation_error",
                    "error_detail": validation_message,
                }
            )
            continue

        source_row = input_rows.get(custom_id, {})
        parsed_records.append(parsed)

        country_evidence = parsed.get("country_evidence", [])
        if country_evidence:
            for item in country_evidence:
                flattened_rows.append(
                    {
                        "record_key": custom_id,
                        "Title": source_row.get("input", {}).get("Title", ""),
                        "article_url": source_row.get("input", {}).get("article_url", ""),
                        "country_standard": item.get("country", ""),
                        "evidence_span": item.get("evidence_span", ""),
                        "reasoning_label": parsed.get("reasoning_label", ""),
                        "confidence": parsed.get("confidence", ""),
                        "llm_source": args.llm_source,
                    }
                )
        else:
            flattened_rows.append(
                {
                    "record_key": custom_id,
                    "Title": source_row.get("input", {}).get("Title", ""),
                    "article_url": source_row.get("input", {}).get("article_url", ""),
                    "country_standard": "",
                    "evidence_span": "",
                    "reasoning_label": parsed.get("reasoning_label", ""),
                    "confidence": parsed.get("confidence", ""),
                    "llm_source": args.llm_source,
                }
            )

    parsed_jsonl_path = Path(args.parsed_jsonl)
    candidates_csv_path = Path(args.candidates_csv)
    invalid_csv_path = Path(args.invalid_csv)

    write_jsonl(parsed_jsonl_path, parsed_records)

    with candidates_csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "record_key",
                "Title",
                "article_url",
                "country_standard",
                "evidence_span",
                "reasoning_label",
                "confidence",
                "llm_source",
            ],
        )
        writer.writeheader()
        writer.writerows(flattened_rows)

    with invalid_csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["record_key", "error_type", "error_detail"],
        )
        writer.writeheader()
        writer.writerows(invalid_rows)

    message(f"Parsed valid batch records: {len(parsed_records)}")
    message(f"Wrote candidate evidence CSV: {candidates_csv_path}")
    message(f"Wrote invalid output log: {invalid_csv_path}")


def build_parser() -> argparse.ArgumentParser:
    batch_dir = ensure_batch_dir()
    parser = argparse.ArgumentParser(
        description="OpenAI Batch runner for unresolved WHO DON country extraction."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="Build OpenAI Batch request JSONL.")
    prepare.add_argument("--input-jsonl", default=str(INPUT_JSONL))
    prepare.add_argument(
        "--requests-jsonl",
        default=str(batch_dir / "who_don_openai_batch_requests.jsonl"),
    )
    prepare.add_argument("--model", default=DEFAULT_MODEL)
    prepare.add_argument(
        "--workflow",
        default="extract",
        choices=["extract", "adjudicate"],
    )
    prepare.add_argument(
        "--endpoint",
        default=DEFAULT_ENDPOINT,
        choices=["/v1/responses", "/v1/chat/completions"],
    )
    prepare.add_argument(
        "--reasoning-effort",
        default="low",
        choices=["none", "low", "medium", "high", "xhigh"],
    )
    prepare.add_argument("--limit", type=int, default=None)
    prepare.add_argument(
        "--record-keys",
        default=None,
        help="Comma-separated record_key values to include in the prepared batch.",
    )
    prepare.set_defaults(func=prepare_requests)

    submit = subparsers.add_parser("submit", help="Upload requests and create a batch.")
    submit.add_argument(
        "--requests-jsonl",
        default=str(batch_dir / "who_don_openai_batch_requests.jsonl"),
    )
    submit.add_argument("--completion-window", default=DEFAULT_COMPLETION_WINDOW)
    submit.add_argument("--model", default=DEFAULT_MODEL)
    submit.add_argument(
        "--workflow",
        default="extract",
        choices=["extract", "adjudicate"],
    )
    submit.add_argument(
        "--endpoint",
        default=DEFAULT_ENDPOINT,
        choices=["/v1/responses", "/v1/chat/completions"],
    )
    submit.add_argument(
        "--upload-json",
        default=str(batch_dir / "who_don_openai_batch_upload.json"),
    )
    submit.add_argument(
        "--batch-json",
        default=str(batch_dir / "who_don_openai_batch_job.json"),
    )
    submit.set_defaults(func=submit_batch)

    status = subparsers.add_parser("status", help="Fetch current batch job status.")
    status.add_argument("--batch-id", required=True)
    status.add_argument(
        "--batch-json",
        default=str(batch_dir / "who_don_openai_batch_job.json"),
    )
    status.set_defaults(func=fetch_batch_status)

    download = subparsers.add_parser("download", help="Download completed batch files.")
    download.add_argument(
        "--batch-json",
        default=str(batch_dir / "who_don_openai_batch_job.json"),
    )
    download.add_argument(
        "--output-jsonl",
        default=str(batch_dir / "who_don_openai_batch_output.jsonl"),
    )
    download.add_argument(
        "--error-jsonl",
        default=str(batch_dir / "who_don_openai_batch_error.jsonl"),
    )
    download.set_defaults(func=download_batch_outputs)

    parse = subparsers.add_parser("parse", help="Parse and validate batch outputs.")
    parse.add_argument("--input-jsonl", default=str(INPUT_JSONL))
    parse.add_argument(
        "--output-jsonl",
        default=str(batch_dir / "who_don_openai_batch_output.jsonl"),
    )
    parse.add_argument(
        "--parsed-jsonl",
        default=str(batch_dir / "who_don_openai_batch_parsed.jsonl"),
    )
    parse.add_argument(
        "--candidates-csv",
        default=str(batch_dir / "who_don_openai_country_candidates.csv"),
    )
    parse.add_argument(
        "--invalid-csv",
        default=str(batch_dir / "who_don_openai_invalid_outputs.csv"),
    )
    parse.add_argument("--llm-source", default="openai_batch")
    parse.set_defaults(func=parse_batch_outputs)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        args.func(args)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
