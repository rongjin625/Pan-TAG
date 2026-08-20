"""Shared configurable GFF/GTF attribute handling for Pan-TAG."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class AnnotationSchema:
    mode: str
    gene_feature_types: tuple[str, ...]
    pseudogene_feature_types: tuple[str, ...]
    transcript_feature_types: tuple[str, ...]
    gene_id_keys: tuple[str, ...]
    gene_name_keys: tuple[str, ...]
    transcript_id_keys: tuple[str, ...]
    transcript_parent_gene_keys: tuple[str, ...]
    biotype_keys: tuple[str, ...]
    coding_filter_mode: str
    coding_values: tuple[str, ...]
    excluded_attribute_terms: tuple[str, ...]
    pseudogene_values: tuple[str, ...]
    strip_id_prefixes: tuple[str, ...]
    fallback_gene_name_to_id: bool


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            key, value = item, ""
        attributes[key.strip()] = value.strip().strip('"').strip("'")
    return attributes


def first_value(attributes: dict[str, str], keys: tuple[str, ...]) -> str:
    for key in keys:
        value = attributes.get(key, "").strip()
        if value:
            return value.split(",", 1)[0].strip()
    return ""


def all_values(
    attributes: dict[str, str],
    keys: tuple[str, ...],
    prefixes: tuple[str, ...],
) -> tuple[str, ...]:
    """Return all unique, normalized values for the configured keys."""
    values: list[str] = []
    seen: set[str] = set()
    for key in keys:
        for raw_value in attributes.get(key, "").split(","):
            value = normalize_id(raw_value, prefixes)
            if value and value not in seen:
                seen.add(value)
                values.append(value)
    return tuple(values)


def normalize_id(value: str, prefixes: tuple[str, ...]) -> str:
    result = value.strip()
    changed = True
    while changed:
        changed = False
        for prefix in prefixes:
            if prefix and result.startswith(prefix):
                result = result[len(prefix):]
                changed = True
                break
    return result


def load_schema(config_path: str | Path) -> AnnotationSchema:
    with Path(config_path).open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    annotation = config["annotation"]

    def values(key: str) -> tuple[str, ...]:
        configured = annotation[key]
        if isinstance(configured, str):
            configured = [configured]
        result = tuple(str(value) for value in configured if str(value))
        if not result:
            raise ValueError(f"annotation.{key} must contain at least one value")
        return result

    coding_filter_mode = str(
        annotation["coding_filter_mode"]
    ).lower()
    if coding_filter_mode not in {"auto", "biotype", "translated_protein"}:
        raise ValueError(
            "annotation.coding_filter_mode must be auto, biotype or translated_protein"
        )

    return AnnotationSchema(
        mode=str(annotation["mode"]).lower(),
        gene_feature_types=values("gene_feature_types"),
        pseudogene_feature_types=values("pseudogene_feature_types"),
        transcript_feature_types=values("transcript_feature_types"),
        gene_id_keys=values("gene_id_keys"),
        gene_name_keys=values("gene_name_keys"),
        transcript_id_keys=values("transcript_id_keys"),
        transcript_parent_gene_keys=values("transcript_parent_gene_keys"),
        biotype_keys=values("biotype_keys"),
        coding_filter_mode=coding_filter_mode,
        coding_values=values("coding_values"),
        excluded_attribute_terms=values("excluded_attribute_terms"),
        pseudogene_values=values("pseudogene_values"),
        strip_id_prefixes=values("strip_id_prefixes"),
        fallback_gene_name_to_id=bool(annotation["fallback_gene_name_to_id"]),
    )


def gene_id(attributes: dict[str, str], schema: AnnotationSchema) -> str:
    return normalize_id(first_value(attributes, schema.gene_id_keys), schema.strip_id_prefixes)


def gene_name(attributes: dict[str, str], schema: AnnotationSchema, identifier: str) -> tuple[str, bool]:
    name = first_value(attributes, schema.gene_name_keys)
    if name:
        return name, False
    if schema.fallback_gene_name_to_id:
        return identifier, True
    return "", False


def transcript_id(attributes: dict[str, str], schema: AnnotationSchema) -> str:
    return normalize_id(first_value(attributes, schema.transcript_id_keys), schema.strip_id_prefixes)


def transcript_ids(attributes: dict[str, str], schema: AnnotationSchema) -> tuple[str, ...]:
    return all_values(attributes, schema.transcript_id_keys, schema.strip_id_prefixes)


def transcript_parent(attributes: dict[str, str], schema: AnnotationSchema) -> str:
    return normalize_id(
        first_value(attributes, schema.transcript_parent_gene_keys), schema.strip_id_prefixes
    )


def is_pseudogene(attributes: dict[str, str], schema: AnnotationSchema) -> bool:
    configured = {value.lower() for value in schema.pseudogene_values}
    values = [attributes.get(key, "").lower() for key in schema.biotype_keys]
    return any(value in configured or "pseudogene" in value for value in values if value)
