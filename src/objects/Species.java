package objects;

import core.Utils;

import java.util.ArrayList;
import java.util.List;

public class Species {

    private final String name;
    private final List<Gene> genes;

    public Species(String name) {
        this.name = name;
        this.genes = new ArrayList<>();
    }

    public double findBestMatchDistance(Gene comparedGene) {
        double currentBest = Double.MAX_VALUE;
        for (Gene gene : this.getGenes()) {
            double distance = Utils.getDistanceBetweenGenes(comparedGene, gene);
            if (distance < currentBest) {
                currentBest = distance;
            }
        }
        return currentBest;
    }

    public List<Gene> getGenes() {
        return genes;
    }

    public String getName() {
        return name;
    }
}
