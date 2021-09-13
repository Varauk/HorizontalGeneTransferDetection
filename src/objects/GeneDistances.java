package objects;

import objects.Matrix.MatrixFormat;

import java.util.ArrayList;
import java.util.List;

public class GeneDistances {

    private List<Species> species;
    private List<Gene> genes;
    private Matrix<Double> distances;
    private String treeName;

    /**
     * Constructor declares ArrayList genes and ArrayList distances and sets the attribute treeName.
     *
     * @param treeName Identifier of the root of the gene tree as String
     * @param format   Determines, if the matrix is an upper, lower or symmetrical matrix
     */
    public GeneDistances(MatrixFormat format, String treeName) {
        this.species = new ArrayList<>();
        this.genes = new ArrayList<Gene>();
        this.distances = new Matrix<Double>(format);
        this.treeName = treeName;
    }

    /**
     * Getter.
     *
     * @return ArrayList of Gene
     */
    public List<Gene> getGenes() {
        return this.genes;
    }

    /**
     * Getter.
     *
     * @return ArrayList of ArrayLists of Doubles represents a lower triangle matrix
     */
    public Matrix<Double> getDistances() {
        return this.distances;
    }

    /**
     * Getter.
     *
     * @return String with tree name
     */
    public String getTreeName() {
        return this.treeName;
    }

    public void setSpecies(List<Species> species) {
        this.species = species;
    }

    public List<Species> getSpecies() {
        return species;
    }

    public void setGenes(List<Gene> genes) {
        this.genes = genes;
    }
}
