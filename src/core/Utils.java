package core;

import objects.Gene;
import objects.GeneTuple;
import objects.Species;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

public abstract class Utils {
    public static double getDistanceBetweenGenes(Gene gene, Gene comparedGene) {
        return gene.getDistances().getOrDefault(comparedGene, Double.MAX_VALUE);
    }

    public static List<Species> getOrderedSpecies(Map<Species, Double> map) {
        List<Species> orderedSpecies = new ArrayList<>();
        //map.keySet().length^2 runtime
        for (int i = 0; i < map.keySet().size(); i++) {
            Species next = null;
            Double nextDouble = Double.MAX_VALUE;
            for (Species sForOrder : map.keySet()) {
                //smaller than next int AND not already in the orderedSpecies list
                if (map.get(sForOrder) <= nextDouble && orderedSpecies.stream().filter(o -> o.getName().equals(sForOrder.getName())).findFirst().orElse(null) == null) {
                    nextDouble = map.get(sForOrder);
                    next = sForOrder;
                }
            }
            orderedSpecies.add(next);
        }
        return orderedSpecies;
    }

    public static void evaluateResults(List<Species> allSpeciesAfterRunningAlgorithm, List<Gene> genesWithRCBRelationsExtractedFromInput) {

        Logger log = Logger.getGlobal();
        if (genesWithRCBRelationsExtractedFromInput == null || genesWithRCBRelationsExtractedFromInput.size() == 0) {
            log.warning("Skipping evaluation.");
            return;
        }
        int recognizedHGTs = 0;
        int missedHGTs = 0;
        int wrongHGTs = 0;
        int bothEqual = 0;
        int correctHigher = 0;
        int incorrectHigher = 0;
        int correctLesser = 0;
        int incorrectLesser = 0;
        int lesserAsEqual = 0;
        int higherAsEqual = 0;
        for (Species s : allSpeciesAfterRunningAlgorithm) {
            for (Gene gene : s.getGenes()) {
                Gene matchingGeneFromInput = genesWithRCBRelationsExtractedFromInput.stream().filter(g -> g.getSpeciesAndGeneIdentifier().equals(gene.getSpeciesAndGeneIdentifier())).findFirst().orElse(null);
                if (matchingGeneFromInput != null) {
                    List<GeneTuple> realRelations = matchingGeneFromInput.getRelationsToOtherGenesInItsTree();
                    if (realRelations.size() != gene.getRelationsToOtherGenesInItsTree().size()) {
                        log.warning("Relation list size mismatch, skipping " + gene.getSpeciesAndGeneIdentifier());
                        continue;
                    }

                    for (int i = 0; i < gene.getRelationsToOtherGenesInItsTree().size(); i++) {
                        GeneTuple foundRelation = gene.getRelationsToOtherGenesInItsTree().get(i);
                        GeneTuple realRelation = null;
                        for (GeneTuple potentialRealRelation : realRelations) {
                            if (potentialRealRelation.containsGenes(foundRelation.first.getSpeciesAndGeneIdentifier(), foundRelation.second.getSpeciesAndGeneIdentifier())) {
                                realRelation = potentialRealRelation;
                            }
                        }

                        if (realRelation == null) {
                            log.severe("Didnt find a matching real relation!");
                            return;
                        }

                        if (foundRelation.type == GeneTuple.RelationType.EQUAL) {
                            //algorithm didn't find anything, which is to some degree expected
                            if (realRelation.type != GeneTuple.RelationType.EQUAL) {
                                missedHGTs++;
                                if (realRelation.type == GeneTuple.RelationType.LESSER) {
                                    lesserAsEqual++;
                                } else if (realRelation.type == GeneTuple.RelationType.HIGHER) {
                                    higherAsEqual++;
                                }
                            } else {
                                bothEqual++;
                            }

                        } else if (foundRelation.type == GeneTuple.RelationType.LESSER) {
                            if (realRelation.type == GeneTuple.RelationType.LESSER) {
                                correctLesser++;
                                recognizedHGTs++;
                            } else {
                                wrongHGTs++;
                                incorrectLesser++;
                                if (realRelation.type == GeneTuple.RelationType.HIGHER) {
                                    //incorrectLesser++;
                                }
                                //log.warning("Wrong LESSER HGT found with gene " + realRelation.first.getSpeciesAndGeneIdentifier() + " and " + realRelation.second.getSpeciesAndGeneIdentifier() + ", gene tree " + gene.getGeneTree());
                            }

                        } else if (foundRelation.type == GeneTuple.RelationType.HIGHER) {
                            if (realRelation.type == GeneTuple.RelationType.HIGHER) {
                                correctHigher++;
                                recognizedHGTs++;
                            } else {
                                wrongHGTs++;
                                incorrectHigher++;
                                if (realRelation.type == GeneTuple.RelationType.LESSER) {
                                    //incorrectHigher++;
                                }
                                //log.warning("Wrong HIGHER HGT found with gene " + realRelation.first.getSpeciesAndGeneIdentifier() + " and " + realRelation.second.getSpeciesAndGeneIdentifier() + ", gene tree " + gene.getGeneTree());
                            }
                        }
                    }
                } else {
                    log.warning("Couldn't find a matching gene for " + gene.getSpeciesAndGeneIdentifier() + " in the rcb relation list extracted from the input file.");
                }
            }
        }
        correctHigher = correctHigher / 2;
        incorrectHigher = incorrectHigher / 2;
        correctLesser = correctLesser / 2;
        incorrectLesser = incorrectLesser / 2;
        lesserAsEqual = lesserAsEqual / 2;
        higherAsEqual = higherAsEqual / 2;
        recognizedHGTs = recognizedHGTs / 2;
        wrongHGTs = wrongHGTs / 2;
        int totalFound = recognizedHGTs + wrongHGTs;

        //**START true amount of HGTs in RELATIONS block
        int totalLesser = 0;
        int totalHigher = 0;

        for (Gene gene : genesWithRCBRelationsExtractedFromInput) {
            for (GeneTuple relation : gene.getRelationsToOtherGenesInItsTree()) {
                if (relation.type == GeneTuple.RelationType.LESSER) {
                    totalLesser++;
                } else if (relation.type == GeneTuple.RelationType.HIGHER) {
                    totalHigher++;
                }
            }
        }

        totalLesser = totalLesser / 2;
        totalHigher = totalHigher / 2;
        //**END

        //**START all distances between all genes of two species
        /*final String species1 = "SE007";
        final String species2 = "SE008";
        int[] bestMatchPositions = new int[allSpeciesAfterRunningAlgorithm.size()];
        log.warning("Distances between genes of species " + species1 + " and " + species2);
        List<String> stringList = new ArrayList<>();
        for (Species s : allSpeciesAfterRunningAlgorithm) {
            if (s.getName().equals(species1)) {
                for (Gene gene : s.getGenes()) {
                    for (Gene partnerGene : gene.getDistances().keySet()) {
                        if (partnerGene.getParent().getName().equals(species2)) {
                            //Gene trueGene = genesWithRCBRelationsExtractedFromInput.stream().filter(g1 -> g1.getSpeciesAndGeneIdentifier().equals(gene.getSpeciesAndGeneIdentifier())).findFirst().orElse(null);
                            //GeneTuple.RelationType type = trueGene.getRelationsToOtherGenesInItsTree().stream().filter(tuple -> tuple.containsGenes(gene, partnerGene)).findFirst().orElse(null).type;
                            GeneTuple.RelationType type = genesWithRCBRelationsExtractedFromInput.stream().filter(g1 -> g1.getSpeciesAndGeneIdentifier().equals(gene.getSpeciesAndGeneIdentifier())).findFirst().orElse(null).getRelationsToOtherGenesInItsTree().stream().filter(tuple -> tuple.containsGenes(gene, partnerGene)).findFirst().orElse(null).type;
                            stringList.add(gene.getDistances().get(partnerGene) + " (" + gene.getSpeciesAndGeneIdentifier() + " with " + partnerGene.getSpeciesAndGeneIdentifier() + ", " + type + ", " + gene.getGeneTree() + ")");
                            //stringList.add("$" + gene.getDistances().get(partnerGene) + "$ & " + gene.getGeneTree() + " & " + type + "\\\\");
                        }
                    }
                    List<DistanceSpecieslistTuple> list = new ArrayList<>();
                    DistanceSpecieslistTuple.addDistanceSpecieslistTuplesToList(list, gene.getBestMatchDistances());
                    for (int i = 0; i < list.size(); i++) {
                        if (list.get(i).species.stream().filter(partnerSpecies -> partnerSpecies.getName().equals(species2)).findFirst().orElse(null) != null) {
                            bestMatchPositions[i]++;
                            break;
                        }
                    }
                }
            }
        }

        stringList.sort(String::compareTo);
        for (String s : stringList) {
            log.warning(s);
            System.out.println(s);
        }*/
        //**END


        //divided by 2 because both directions are counted here
        DecimalFormat df = new DecimalFormat("###.##");
        log.warning("Total found: $" + totalFound + "$. Both Equal: " + (double) bothEqual / 2);
        log.warning("Correctly recognized HGTs: $" + (double) recognizedHGTs + "$ $(" + df.format(((double) recognizedHGTs) / totalFound * 100) + "\\%)$");
        log.warning("Wrong HGTs (false positives): $" + (double) wrongHGTs + "$ $(" + df.format(((double) wrongHGTs) / totalFound * 100) + "\\%)$");
        log.warning("Missed HGTs (false negatives): $" + (double) missedHGTs / 2 + "$");
        log.warning("Total lesser: " + totalLesser + ", Total higher: " + totalHigher + ", Total existing: " + (totalLesser + totalHigher));
        //System.out.println("& $" + totalFound + "$ & $" + (double) recognizedHGTs + "$ $(" + df.format(((double) recognizedHGTs) / totalFound * 100) + "\\%)$ & $" + (double) wrongHGTs + "$ $(" + df.format(((double) wrongHGTs) / totalFound * 100) + "\\%)$ & $" + (double) missedHGTs / 2 + "$ \\\\");
        System.out.println("& $" + totalFound + "$ & $" + recognizedHGTs + "$ $(" + df.format(((double) recognizedHGTs) / totalFound * 100) + "\\%)$ & $" + wrongHGTs + "$ $(" + df.format(((double) wrongHGTs) / totalFound * 100) + "\\%)$ & $" + (totalLesser + totalHigher - recognizedHGTs) + "$ \\\\");
    }
}
