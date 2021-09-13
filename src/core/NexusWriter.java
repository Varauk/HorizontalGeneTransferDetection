package core;

import objects.Gene;
import objects.GeneDistances;
import objects.GeneTuple;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Logger;

public class NexusWriter {
    public static void writeHgtsToFile(String filename, HashMap<String, GeneDistances> map, List<GeneTuple> hgts) {
        Logger log = Logger.getGlobal();
        log.warning("Writing to file " + filename + "...");
        try {
            FileWriter file = new FileWriter(filename, false);
            file.append("#NEXUS\nBEGIN HGTRELATIONS;\n[matrices each describing the type of horizontal gene transfer between genes]\n");
            file.append("\t[LESSER=-1, EQUALORDEFAULT=0, GREATER=1]\n");

            for (String genetree : map.keySet()) {
                file.append("\thgtrelation\n\t\tname=" + genetree + " triangle=both\n");
                List<Gene> genes = map.get(genetree).getGenes();
                //genes.sort(Comparator.comparing(Gene::getSpeciesAndGeneIdentifier)); //looks better, but without this the order is the same as in the input nexus file
                for (Gene gene : genes) {
                    file.append("\t\t\t" + gene.getSpeciesAndGeneIdentifier() + " ");
                    for (Gene comparedGene : genes) {
                        if (gene.equals(comparedGene)) {
                            file.append(" 0 ");
                            gene.getRelationsToOtherGenesInItsTree().add(new GeneTuple(gene, gene, GeneTuple.RelationType.EQUAL));
                            continue;
                        }
                        GeneTuple foundHGT = hgts.stream().filter(tuple -> tuple.containsGenes(gene, comparedGene)).findFirst().orElse(null);
                        if (foundHGT == null) {
                            file.append(" 0 ");
                            gene.getRelationsToOtherGenesInItsTree().add(new GeneTuple(gene, comparedGene, GeneTuple.RelationType.EQUAL));
                        } else {
                            if (foundHGT.type == GeneTuple.RelationType.HIGHER) {
                                file.append(" 1 ");
                                gene.getRelationsToOtherGenesInItsTree().add(new GeneTuple(gene, comparedGene, GeneTuple.RelationType.HIGHER));
                            } else if (foundHGT.type == GeneTuple.RelationType.LESSER) {
                                file.append("-1 ");
                                gene.getRelationsToOtherGenesInItsTree().add(new GeneTuple(gene, comparedGene, GeneTuple.RelationType.LESSER));
                            } else {
                                file.append(" 0 ");
                                gene.getRelationsToOtherGenesInItsTree().add(new GeneTuple(gene, comparedGene, GeneTuple.RelationType.EQUAL));
                            }
                        }
                    }
                    file.append("\n");
                }
                file.append("\t\t\t;\n");
            }
            file.append("END;");
            file.close();
        } catch (IOException e) {
            log.severe("Couldn't write to file " + filename);
            e.printStackTrace();
        }
        log.warning("Export completed!");
    }
}
