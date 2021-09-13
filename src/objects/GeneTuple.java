package objects;

public class GeneTuple {
    public final Gene first;
    public final Gene second;
    public final RelationType type;

    public GeneTuple(Gene first, Gene second, RelationType type) {
        if (first.getSpeciesAndGeneIdentifier().compareTo(second.getSpeciesAndGeneIdentifier()) < 0) {
            this.first = first;
            this.second = second;
        } else {
            this.first = second;
            this.second = first;
        }
        this.type = type;
    }

    public GeneTuple(Gene first, Gene second) {
        this(first, second, RelationType.EQUAL);
    }

    public boolean containsGenes(Gene geneToCheck1, Gene geneToCheck2) {
        return (first.getSpeciesAndGeneIdentifier().equals(geneToCheck1.getSpeciesAndGeneIdentifier()) && second.getSpeciesAndGeneIdentifier().equals(geneToCheck2.getSpeciesAndGeneIdentifier()))
                || (first.getSpeciesAndGeneIdentifier().equals(geneToCheck2.getSpeciesAndGeneIdentifier()) && second.getSpeciesAndGeneIdentifier().equals(geneToCheck1.getSpeciesAndGeneIdentifier()));
    }

    public boolean containsGenes(String geneToCheck1, String geneToCheck2) {
        return (first.getSpeciesAndGeneIdentifier().equals(geneToCheck1) && second.getSpeciesAndGeneIdentifier().equals(geneToCheck2)) || (first.getSpeciesAndGeneIdentifier().equals(geneToCheck2) && second.getSpeciesAndGeneIdentifier().equals(geneToCheck1));
    }

    public enum RelationType {
        HIGHER, EQUAL, LESSER
    }
}
