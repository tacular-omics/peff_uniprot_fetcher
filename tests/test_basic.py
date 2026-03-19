"""Basic smoke test for the package."""

import peff_uniprot_fetcher


def test_package_has_public_api():
    assert hasattr(peff_uniprot_fetcher, "fetch_peff")
    assert hasattr(peff_uniprot_fetcher, "fetch_peff_to_file")
    assert hasattr(peff_uniprot_fetcher, "write_peff")
