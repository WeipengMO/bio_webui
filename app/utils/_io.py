def parse_gene_input(genes, remove_duplicates=True):
    genes = genes.split()
    seen = set()
    if remove_duplicates:
        genes = [
            gene.strip() for gene in genes 
                if gene.strip() not in seen and not seen.add(gene.strip())]  # Remove duplicates
    else:
        genes = [gene.strip() for gene in genes]
        
    return genes