package objects;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class DistanceSpecieslistTuple {
    public final double distance;
    public final List<Species> species;

    /*public DistanceGenelistTuple(double distance, List<Gene> genes) {
        this.distance = distance;
        this.genes = genes;
    }*/

    public DistanceSpecieslistTuple(double distance, Species firstSpecies) {
        this.distance = distance;
        List<Species> species = new ArrayList<>();
        species.add(firstSpecies);
        this.species = species;
    }

    /**
     * Adds every species in the map to the list, and groups them in DistanceSpecieslistTuples by distance
     *
     * @param list List of DistanceSpecieslistTuple which we want to add stuff to
     * @param map  Map of SpeciesXDouble, usually coming from Gene.getBestMatchDistances()
     */
    public static void addDistanceSpecieslistTuplesToList(List<DistanceSpecieslistTuple> list, Map<Species, Double> map) {
        for (Species species : map.keySet()) {
            DistanceSpecieslistTuple existingTuple = list.stream().filter(tuple -> tuple.distance == map.get(species)).findFirst().orElse(null);
            if (existingTuple == null) {
                list.add(new DistanceSpecieslistTuple(map.get(species), species));
            } else {
                existingTuple.species.add(species);
            }
        }
        list.sort(Comparator.comparing(tuple -> tuple.distance));
    }
}
