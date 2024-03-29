"""Utilities to parse the PSI-MI TAB file format."""

from typing import Container, Optional


def parse(entry: str) -> dict[str, tuple[tuple[str, Optional[str]], ...]]:
    """
    Parse an entry in PSI-MI TAB format.

    Args:
        entry: A PSI-MI TAB formatted entry, corresponding to
            'namespace:"identifier"(term)'.

    Returns:
        A keyed representation of the entry, separating namespace, identifiers
        and terms
    """
    # Parse components PSI-MI TAB terms.
    if entry == "-":
        return {}
    else:
        values: dict[str, list[tuple[str, Optional[str]]]] = {}
        for namespace, identifier_term in (namespace_identifier.split(
                ":", 1) for namespace_identifier in entry.split("|")):

            if "(" in identifier_term and ")" in identifier_term:
                identifier = identifier_term[:identifier_term.find("(")].strip(
                    '"')

                term = identifier_term[identifier_term.find("(") +
                                       1:identifier_term.find(")")].strip('"')
            else:
                identifier = identifier_term.strip('"')
                term = None

            if namespace in values:
                values[namespace].append((identifier, term))
            else:
                values[namespace] = [(identifier, term)]

        return {
            ns: tuple(identifiers_terms)
            for ns, identifiers_terms in values.items()
        }


def get_identifiers_from_namespace(entry: str,
                                   namespace: str) -> tuple[str, ...]:
    """
    Returns identifiers from a namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.

    Returns:
        Identifiers from a namespace.
    """
    # Extract PSI-MI TAB identifiers corresponding to a specific namespace.
    return tuple(identifier
                 for identifier, _ in parse(entry).get(namespace, [("", None)])
                 if identifier)


def get_terms_from_namespace(entry: str, namespace: str) -> tuple[str, ...]:
    """
    Returns terms from a namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.

    Returns:
        Terms from a namespace.
    """
    # Extract PSI-MI TAB terms corresponding to a specific namespace.
    return tuple(
        term for _, term in parse(entry).get(namespace, [("", None)]) if term)


def namespace_has_identifier(entry: str, namespace: str,
                             identifier: int | str) -> bool:
    """
    Returns whether an identifier exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        identifier: The identifier to search.

    Returns:
        True, if identifier is in namespace, else False.
    """
    # Indicate whether a namespace has a specific PSI-MI TAB identifier.
    return str(identifier) in get_identifiers_from_namespace(entry, namespace)


def namespace_has_term(entry: str, namespace: str, term: int | str) -> bool:
    """
    Returns whether a term exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        term: The term to search.

    Returns:
        True, if term is in namespace, else False.
    """
    # Indicate whether a namespace has a specific PSI-MI TAB term.
    return str(term) in get_terms_from_namespace(entry, namespace)


def namespace_has_any_identifier_from(
        entry: str, namespace: str, identifiers: Container[int | str]) -> bool:
    """
    Returns whether any identifier exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        identifiers: The identifiers to search.

    Returns:
        True, if any identifier is in namespace, else False.
    """
    # Indicate whether a namespace has any from specific PSI-MI TAB identifiers.
    return any(
        identifier in identifiers
        for identifier in get_identifiers_from_namespace(entry, namespace))


def namespace_has_any_term_from(entry: str, namespace: str,
                                terms: Container[int | str]) -> bool:
    """
    Returns whether any term exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        terms: The terms to search.

    Returns:
        True, if any term is in namespace, else False.
    """
    # Indicate whether a namespace has any from specific PSI-MI TAB terms.
    return any(
        term in terms for term in get_terms_from_namespace(entry, namespace))
