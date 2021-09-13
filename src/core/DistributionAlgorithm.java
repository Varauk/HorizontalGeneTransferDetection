package core;

import objects.Counter;
import objects.DistanceSpecieslistTuple;
import objects.Gene;
import objects.GeneTuple;
import objects.Species;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

public abstract class DistributionAlgorithm {

    public static List<GeneTuple> runAlgorithm(List<Species> allSpecies, Double percentage, boolean multithreading, GeneTuple.RelationType runMode) {
        Logger log = Logger.getGlobal();
        log.info("Algorithm started");
        List<GeneTuple> allFoundHGTs = new ArrayList<>();
        int countGenes = 0;
        for (Species s : allSpecies) {
            countGenes += s.getGenes().size();
        }
        final int countedGenes = countGenes;

        final Counter counter = new Counter();

        Thread status = new Thread(() -> {
            DecimalFormat df = new DecimalFormat("###.##");
            while (true) {
                try {
                    TimeUnit.SECONDS.sleep(30);
                    if (counter.value < countedGenes) {
                        log.warning("Processed " + counter.value + "/" + countedGenes + " genes (" + df.format(((double) counter.value / countedGenes) * 100) + "%)");
                    } else {
                        return;
                    }
                } catch (InterruptedException e) {
                    return;
                }
            }
        });
        status.start();

        //***start main computing part
        if (multithreading) {
            CountDownLatch latch = new CountDownLatch(allSpecies.size());
            for (Species s : allSpecies) {
                Thread thread = new Thread(() -> {
                    evaluateSpecies(s, allSpecies, percentage, counter, runMode);
                    latch.countDown();
                });
                thread.start();
            }

            //wait for all threads to finish
            try {
                latch.await();
            } catch (InterruptedException e) {
                log.severe("Multithreading encountered an error, try to disable it with -disableMultithreading");
                System.exit(-1);
                //e.printStackTrace();
            }
        } else {
            for (Species s : allSpecies) {
                evaluateSpecies(s, allSpecies, percentage, counter, runMode);
            }
        }
        //***end main computing part

        status.interrupt();

        //***start: evaluate every gene.getPotentialCandidate... list and add all HGTs to a list
        for (Species s : allSpecies) {
            for (Gene g : s.getGenes()) {
                for (Gene potentialCandidate : g.getPotentialCandidateGenesLesser()) {
                    if (potentialCandidate.getPotentialCandidateGenesLesser().stream().filter(testGene -> testGene.equals(g)).findFirst().orElse(null) != null) {
                        allFoundHGTs.add(new GeneTuple(g, potentialCandidate, GeneTuple.RelationType.LESSER));
                    }
                    potentialCandidate.getPotentialCandidateGenesLesser().remove(g);
                }

                for (Gene potentialCandidate : g.getPotentialCandidateGenesHigher()) {
                    if (potentialCandidate.getPotentialCandidateGenesHigher().stream().filter(testGene -> testGene.equals(g)).findFirst().orElse(null) != null) {
                        allFoundHGTs.add(new GeneTuple(g, potentialCandidate, GeneTuple.RelationType.HIGHER));
                    }
                    potentialCandidate.getPotentialCandidateGenesHigher().remove(g);
                }
            }
        }
        //***end: create hgt list


        if (allFoundHGTs.size() > 0) {
            allFoundHGTs.sort((t1, t2) -> {
                if (t1.first.getGeneTree().compareTo(t2.first.getGeneTree()) < 0) {
                    return -1;
                } else if (t1.first.getGeneTree().compareTo(t2.first.getGeneTree()) > 0) {
                    return 1;
                } else {
                    if (t1.first.getSpeciesAndGeneIdentifier().compareTo(t2.first.getSpeciesAndGeneIdentifier()) < 0) {
                        return -1;
                    } else if (t1.first.getSpeciesAndGeneIdentifier().compareTo(t2.first.getSpeciesAndGeneIdentifier()) > 0) {
                        return 1;
                    } else {
                        if (t1.second.getSpeciesAndGeneIdentifier().compareTo(t2.second.getSpeciesAndGeneIdentifier()) < 0) {
                            return -1;
                        } else {
                            return 1;
                        }
                    }
                }
            });

            log.info("Printing results:");
            for (GeneTuple tuple : allFoundHGTs) {
                if (!tuple.first.getGeneTree().equals(tuple.second.getGeneTree())) {
                    log.severe("Unexpected gene tree mismatch");
                }
                log.info(tuple.first.getSpeciesAndGeneIdentifier() + " with " + tuple.second.getSpeciesAndGeneIdentifier() + ". Type: " + tuple.type + ". Gene tree: " + tuple.first.getGeneTree() + ".");
            }
        }
        log.warning("Algorithm finished, found " + allFoundHGTs.size() + " matches.");
        return allFoundHGTs;
    }

    /**
     * Evaluate a single species and all its genes by finding potential HGT candidates
     *
     * @param s          species
     * @param allSpecies List of all species (iterating over this)
     * @param percentage percentage to look for in the species distribution, for example smallest 5%
     * @param counter    Counter used for the status update thread
     */
    private static void evaluateSpecies(Species s, List<Species> allSpecies, double percentage, Counter counter, GeneTuple.RelationType runMode) {
        Logger log = Logger.getGlobal();
        //***start: preparation
        final int minimumAmountOfGenes = 2;
        if (s.getGenes().size() < minimumAmountOfGenes) {
            log.warning("Skipping species " + s.getName() + " because it has less than " + minimumAmountOfGenes + " genes.");
            return;
        }

        //save the bestMatchDistances map for every gene of species s
        for (Gene gene : s.getGenes()) {
            Map<Species, Double> bestMatchDistances = new HashMap<>();
            for (Species sForBestMatch : allSpecies) {
                if (sForBestMatch != s) {
                    bestMatchDistances.put(sForBestMatch, sForBestMatch.findBestMatchDistance(gene));
                }
            }
            gene.setBestMatchDistances(bestMatchDistances);
        }
        //***end: preparation

        //***start: main part, compare each gene to other genes
        for (Gene gene : s.getGenes()) {
            //needed for distribution computing only. Already here for performance reasons
            List<Species> orderedSpeciesOfCurrentGene = Utils.getOrderedSpecies(gene.getBestMatchDistances()); //this was once done using the updated genemap using the distance of the looked at candidateGene.
            for (Species potentialCandidateSpecies : gene.getBestMatchDistances().keySet()) {

                //***start: distribution of potentialCandidateSpecies
                int[] distribution = new int[allSpecies.size()]; //allSpecies.size is an upper limit
                int distributionTotalValueCount = 0;
                for (Gene differentGeneOfRootSpecies : s.getGenes()) {

                    List<DistanceSpecieslistTuple> orderedSpeciesOfDifferentGeneGrouped = new ArrayList<>();
                    DistanceSpecieslistTuple.addDistanceSpecieslistTuplesToList(orderedSpeciesOfDifferentGeneGrouped, differentGeneOfRootSpecies.getBestMatchDistances());

                    for (DistanceSpecieslistTuple tuple : orderedSpeciesOfDifferentGeneGrouped) {
                        for (Species potentiallyRemoved : tuple.species) {
                            if (orderedSpeciesOfCurrentGene.
                                    stream()
                                    .filter(species ->
                                            species.getName()
                                                    .equals(potentiallyRemoved
                                                            .getName()))
                                    .findFirst().orElse(null) == null) {
                                tuple.species.remove(potentiallyRemoved);
                            }
                        }
                    }
                    orderedSpeciesOfDifferentGeneGrouped.removeIf(tuple -> tuple.species.size() == 0);

                    for (int i = 0; i < orderedSpeciesOfDifferentGeneGrouped.size(); i++) {
                        if (orderedSpeciesOfDifferentGeneGrouped.get(i).species.stream().filter(species -> species.equals(potentialCandidateSpecies)).findFirst().orElse(null) != null) {
                            distribution[i]++;
                            distributionTotalValueCount++;
                            break;
                        }
                    }
                }

                int computedBorderLesser = 0;
                int howManyElementsDidWeFind = 0;
                for (int i = 0; i < distribution.length; i++) {
                    howManyElementsDidWeFind += distribution[i];
                    if (howManyElementsDidWeFind >= distributionTotalValueCount * percentage) {
                        computedBorderLesser = i;
                        break;
                    }
                }

                int computedBorderHigher = distribution.length;
                howManyElementsDidWeFind = 0;
                for (int i = distribution.length - 1; i >= 0; i--) {
                    howManyElementsDidWeFind += distribution[i];
                    if (howManyElementsDidWeFind >= distributionTotalValueCount * percentage) {
                        computedBorderHigher = i;
                        break;
                    }
                }
                //***end: distribution computing

                //get and sort all genes of potentialCandidateSpecies
                List<Gene> orderedGenesOfPotentialCandidateSpeciesByDistanceToGene = new ArrayList<>(potentialCandidateSpecies.getGenes());
                orderedGenesOfPotentialCandidateSpeciesByDistanceToGene.sort(Comparator.comparing(comparedGene -> Utils.getDistanceBetweenGenes(gene, comparedGene)));

                //***start: compare their position if they are a result of a hgt
                //iterate to find LESSER hgts
                if (runMode != GeneTuple.RelationType.HIGHER) {
                    for (Gene potentialCandidateGene : orderedGenesOfPotentialCandidateSpeciesByDistanceToGene) {
                        if (!gene.getGeneTree().equals(potentialCandidateGene.getGeneTree())) {
                            continue;
                        }

                        int positionInCurrentGene = getPositionOfSpeciesComparedToGene(gene, potentialCandidateGene, potentialCandidateSpecies);

                        if (positionInCurrentGene < computedBorderLesser) {
                            log.info("Potential candidate for LESSER: " + gene.getSpeciesAndGeneIdentifier() + " -> " + potentialCandidateGene.getSpeciesAndGeneIdentifier());
                            gene.getPotentialCandidateGenesLesser().add(potentialCandidateGene);
                        } else {
                            break; //if the current potentialCandidateGene already doesn't fit the rest also won't fit because we ordered them by distance
                        }
                    }
                }


                Collections.reverse(orderedGenesOfPotentialCandidateSpeciesByDistanceToGene); //when searching for HIGHER relation

                //iterate to find HIGHER hgts
                if (runMode != GeneTuple.RelationType.LESSER) {
                    for (Gene potentialCandidateGene : orderedGenesOfPotentialCandidateSpeciesByDistanceToGene) {
                        if (!gene.getGeneTree().equals(potentialCandidateGene.getGeneTree())) {
                            continue;
                        }

                        int positionInCurrentGene = getPositionOfSpeciesComparedToGene(gene, potentialCandidateGene, potentialCandidateSpecies);

                        if (positionInCurrentGene > computedBorderHigher) {
                            log.info("Potential candidate for HIGHER: " + gene.getSpeciesAndGeneIdentifier() + " -> " + potentialCandidateGene.getSpeciesAndGeneIdentifier());
                            gene.getPotentialCandidateGenesHigher().add(potentialCandidateGene);
                        } else {
                            break; //if the current potentialCandidateGene already doesn't fit the rest also won't fit because we ordered them by distance
                        }
                    }
                }
                //***end: comparison
            }
            counter.value++; //increment processed gene counter for status by one
        }
        //***end: main part
    }

    /**
     * Returns the position of potentialCandidateSpecies when looking at the best match distances of gene,
     * but using the updated position of potentialCandidateGene
     *
     * @param gene                      root gene
     * @param potentialCandidateGene    gene that is a potential hgt candidate
     * @param potentialCandidateSpecies species of the potentialCandidateGene
     * @return int representing the position of potentialCandidateSpecies when looking at the best match distances
     */
    private static int getPositionOfSpeciesComparedToGene(Gene gene, Gene potentialCandidateGene, Species potentialCandidateSpecies) {
        //***start: update the position of potentialCandidateSpecies in gene.getBestMatchDistances according to gene's distance to potentialCandidateGene

        //overwrite the bestMatchDistance of gene to potentialCandidateSpecies with the current distance of the potentialCandidateGene
        Map<Species, Double> currentGeneMap = gene.getBestMatchDistances();
        Double originalValue = currentGeneMap.get(potentialCandidateSpecies);
        currentGeneMap.put(potentialCandidateSpecies, Utils.getDistanceBetweenGenes(gene, potentialCandidateGene));

        List<DistanceSpecieslistTuple> orderedSpeciesOfCurrentGeneGrouped = new ArrayList<>();
        DistanceSpecieslistTuple.addDistanceSpecieslistTuplesToList(orderedSpeciesOfCurrentGeneGrouped, currentGeneMap);

        currentGeneMap.put(potentialCandidateSpecies, originalValue); //restore correct gene map
        int positionInCurrentGene = 0;

        for (int i = 0; i < orderedSpeciesOfCurrentGeneGrouped.size(); i++) {
            if (orderedSpeciesOfCurrentGeneGrouped.get(i).species.stream().filter(species -> species.equals(potentialCandidateSpecies)).findFirst().orElse(null) != null) {
                positionInCurrentGene = i;
                break;
            }
        }
        //***end: return the updated position
        return positionInCurrentGene;
    }
}
