package objects;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Gene {

    /**
     * Unique ID in the parent's gene list
     */
    private final String id;
    private final Species parent;
    private Map<Species, Double> bestMatchDistances;
    private final List<Gene> potentialCandidateGenesLesser;
    private final List<Gene> potentialCandidateGenesHigher;
    private final Map<Gene, Double> distances;
    private final String geneTree;
    private final List<GeneTuple> relationsToOtherGenesInItsTree;
    private String[] relationsToOtherGenesInItsTreeAsString;

    /**
     * Constructor. Adds the new gene to the parent's gene list
     *
     * @param id
     * @param parent
     * @param geneTree
     */
    public Gene(String id, Species parent, String geneTree) {
        this.id = id;
        this.parent = parent;
        this.bestMatchDistances = null;
        this.potentialCandidateGenesLesser = new ArrayList<>();
        this.potentialCandidateGenesHigher = new ArrayList<>();
        this.relationsToOtherGenesInItsTree = new ArrayList<>();
        this.distances = new HashMap<>();
        this.geneTree = geneTree;
        parent.getGenes().add(this);
    }

    public String getId() {
        return id;
    }

    public Species getParent() {
        return parent;
    }

    public Map<Species, Double> getBestMatchDistances() {
        return bestMatchDistances;
    }

    public void setBestMatchDistances(Map<Species, Double> bestMatchDistances) {
        this.bestMatchDistances = bestMatchDistances;
    }

    public List<Gene> getPotentialCandidateGenesLesser() {
        return potentialCandidateGenesLesser;
    }

    public List<Gene> getPotentialCandidateGenesHigher() {
        return potentialCandidateGenesHigher;
    }

    public List<GeneTuple> getRelationsToOtherGenesInItsTree() {
        return relationsToOtherGenesInItsTree;
    }

    public void setRelationsToOtherGenesInItsTreeAsString(String[] relationsToOtherGenesInItsTreeAsString) {
        this.relationsToOtherGenesInItsTreeAsString = relationsToOtherGenesInItsTreeAsString;
    }

    public String[] getRelationsToOtherGenesInItsTreeAsString() {
        return relationsToOtherGenesInItsTreeAsString;
    }

    public Map<Gene, Double> getDistances() {
        return distances;
    }

    public String getSpeciesAndGeneIdentifier() {
        return this.getParent().getName() + "/" + this.getId();
    }

    public String getGeneTree() {
        return this.geneTree;
    }
}
