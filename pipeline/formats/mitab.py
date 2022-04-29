"""Utilities to parse the PSI-MI TAB file format."""


def parse(entry: str) -> dict[str, str]:
    """
    Parse an entry in PSI-MI TAB format.

    Args:
        entry: A PSI-MI TAB formatted entry, corresponding to
            'namespace:"identifier"(term)'.

    Returns:
        A keyed representation of the entry, separating namespace, identifiers
        and terms
    """
    if entry == "-":
        return {}
    else:
        values = {}
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

        return values


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
    return tuple(term for _, term in parse(entry).get(namespace, [("", None)])
                 if term)


def namespace_has_identifier(entry: str, namespace: str,
                             identifier: str) -> bool:
    """
    Returns whether an identifier exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        identifier: The identifier to search.

    Returns:
        True, if identifier is in namespace, else False.
    """
    return str(identifier) in get_identifiers_from_namespace(entry, namespace)


def namespace_has_term(entry: str, namespace: str, term: str) -> bool:
    """
    Returns whether a term exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        term: The term to search.

    Returns:
        True, if term is in namespace, else False.
    """
    return str(term) in get_terms_from_namespace(entry, namespace)


def namespace_has_any_identifier_from(entry: str, namespace: str,
                                      identifiers: list[str]) -> bool:
    """
    Returns whether any identifier exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        identifiers: The identifiers to search.

    Returns:
        True, if any identifier is in namespace, else False.
    """
    return any(
        namespace_has_identifier(entry, namespace, str(identifier))
        for identifier in identifiers)


def namespace_has_any_term_from(entry: str, namespace: str,
                                terms: list[str]) -> bool:
    """
    Returns whether any term exists in namespace.

    Args:
        entry: The PSI-MITAB entry to parse.
        namespace: The namespace to consider.
        terms: The terms to search.

    Returns:
        True, if any term is in namespace, else False.
    """
    return any(
        namespace_has_term(entry, namespace, str(term)) for term in terms)
