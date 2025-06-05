def test_term_uniqueness():
    taxonomy = CorrosionTaxonomy()
    all_terms = []
    for category in taxonomy.metal_terms.values():
        all_terms.extend(category)
    assert len(all_terms) == len(set(all_terms)), "Duplicate terms detected"
