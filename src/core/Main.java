package core;

import objects.Gene;
import objects.GeneTuple;
import objects.Matrix;
import objects.Species;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {
    public static void main(String[] args) {
        System.setProperty("java.util.logging.SimpleFormatter.format", "[%1$tT.%1$tL] %5$s %n");
        Logger log = Logger.getGlobal();
        log.setLevel(Level.WARNING);
        long start = System.currentTimeMillis();

        String filename = null;
        String outputfilename = "output.nex";
        double percentage = 0.05;
        boolean multithreading = true;

        GeneTuple.RelationType runMode = null;

        for (int i = 0; i < args.length; i++) {
            try {
                if (args[i].equalsIgnoreCase("-i")) {
                    filename = args[i + 1];
                    i++;
                } else if (args[i].equalsIgnoreCase("-o")) {
                    outputfilename = args[i + 1];
                    i++;
                } else if (args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("-silent")) {
                    log.setLevel(Level.SEVERE);
                } else if (args[i].equalsIgnoreCase("-v") || args[i].equalsIgnoreCase("-verbose")) {
                    log.setLevel(Level.INFO);
                } else if (args[i].equalsIgnoreCase("-p")) {
                    percentage = Double.parseDouble(args[i + 1]);
                    i++;
                } else if (args[i].equalsIgnoreCase("-disableMultithreading") || args[i].equalsIgnoreCase("-dm")) {
                    multithreading = false;
                } else if (args[i].equalsIgnoreCase("-higher")) {
                    runMode = GeneTuple.RelationType.HIGHER;
                } else if (args[i].equalsIgnoreCase("-lesser")) {
                    runMode = GeneTuple.RelationType.LESSER;
                } else if (args[i].equalsIgnoreCase("-h") || args[i].equalsIgnoreCase("-help") || args[i].equalsIgnoreCase("help")) {
                    log.severe("Available parameters\n" +
                            "-i filename\t\tSpecify input file (Nexus format containing an ALLDISTANCES block)\n" +
                            "-o filename\t\tOutput nexus file to export relation matrix to\n" +
                            "-s -silent\t\tRemoves most command line output\n" +
                            "-v -verbose\t\tAdds additional command line output, prints every potential candidate\n" +
                            "-p percentage\t\tPercentage of the distribution of a species to look for HGTs in\n" +
                            "-higher\t\t\tOnly search for HGTs of the HIGHER relation\n" +
                            "-lesser\t\t\tOnly search for HGTs of the LESSER relation\n" +
                            "-disableMultithreading\tReduces CPU load but also drastically increases runtime\n" +
                            "-h -help\t\tDisplays this help");
                } else {
                    log.warning("Didn't recognize parameter \"" + args[i] + "\", skipping...");
                }
            } catch (IndexOutOfBoundsException e) {
                log.severe("Malformed parameters, somewhere near \"" + args[args.length - 1] + "\".");
                System.exit(1);
            }
        }

        log.warning("Starting tool...");

        if (filename == null) {
            log.severe("You did not enter an input file. Do so with \"-i filename.nex\"!");
            System.exit(1);
        }

        NexusReaderDistances nr = null;
        try {
            nr = new NexusReaderDistances(filename);
        } catch (FileNotFoundException e) {
            log.severe("Couldn't find file " + filename + ", application is exiting!");
            System.exit(1);
        } catch (Exception e) {
            e.printStackTrace();
        }

        List<Species> allSpecies = new ArrayList<>();

        for (String key : nr.getDistanceMap().keySet()) {
            for (Species s : nr.getDistanceMap().get(key).getSpecies()) {
                Species fromList = allSpecies.stream().filter(species -> species.getName().equals(s.getName())).findFirst().orElse(null);
                if (fromList == null) {
                    allSpecies.add(s);
                } else {
                    fromList.getGenes().addAll(s.getGenes());
                }
            }
            nr.getDistanceMap().get(key).getDistances().changeFormat(Matrix.MatrixFormat.BOTH);

            for (int indexOfGeneThatNeedsDistances = 0; indexOfGeneThatNeedsDistances < nr.getDistanceMap().get(key).getGenes().size(); indexOfGeneThatNeedsDistances++) {
                for (int indexOfDistanceValue = 0; indexOfDistanceValue < nr.getDistanceMap().get(key).getGenes().size(); indexOfDistanceValue++) {
                    nr.getDistanceMap().get(key).getGenes().get(indexOfGeneThatNeedsDistances).getDistances().put(
                            nr.getDistanceMap().get(key).getGenes().get(indexOfDistanceValue),
                            nr.getDistanceMap().get(key).getDistances().get(indexOfGeneThatNeedsDistances).get(indexOfDistanceValue));
                }
            }
        }

        allSpecies.sort(Comparator.comparing(Species::getName));
        allSpecies.forEach(s -> s.getGenes().sort(Comparator.comparing(Gene::getId)));
        log.warning("Running for " + filename);
        List<GeneTuple> foundHgts = DistributionAlgorithm.runAlgorithm(allSpecies, percentage, multithreading, runMode);
        NexusWriter.writeHgtsToFile(outputfilename, nr.getDistanceMap(), foundHgts); //when disabling the following evaluate block uncomment this

        // start of debugging block to load an additional file which contains the true relations, only for testing purposes with ALFSim
//        NexusReaderDistances nrOut = null;
//        final String outFileName = "out1000.nex";
//        try {
//            nrOut = new NexusReaderDistances(".\\src\\resources\\" + outFileName);
//        } catch (FileNotFoundException e) {
//            log.severe(e.toString());
//            log.warning("Couldn't find file " + filename + " in .\\src\\resources\\, trying root now...");
//            try {
//                nr = new NexusReaderDistances(filename);
//            } catch (Exception e2) {
//                log.severe("Couldn't find file " + filename + " in root either, application is exiting!");
//                System.exit(1);
//            }
//        } catch (Exception e) {
//            System.out.println(e);
//        }
//
//        List<Species> allSpeciesOut = new ArrayList<>();
//
//        for (String key : nrOut.getDistanceMap().keySet()) {
//            for (Species s : nrOut.getDistanceMap().get(key).getSpecies()) {
//                Species fromList = allSpeciesOut.stream().filter(species -> species.getName().equals(s.getName())).findFirst().orElse(null);
//                if (fromList == null) {
//                    allSpeciesOut.add(s);
//                } else {
//                    fromList.getGenes().addAll(s.getGenes());
//                }
//            }
//            nrOut.getDistanceMap().get(key).getDistances().changeFormat(Matrix.MatrixFormat.BOTH);
//
//            for (int indexOfGeneThatNeedsDistances = 0; indexOfGeneThatNeedsDistances < nrOut.getDistanceMap().get(key).getGenes().size(); indexOfGeneThatNeedsDistances++) {
//                for (int indexOfDistanceValue = 0; indexOfDistanceValue < nrOut.getDistanceMap().get(key).getGenes().size(); indexOfDistanceValue++) {
//                    nrOut.getDistanceMap().get(key).getGenes().get(indexOfGeneThatNeedsDistances).getDistances().put(
//                            nrOut.getDistanceMap().get(key).getGenes().get(indexOfDistanceValue),
//                            nrOut.getDistanceMap().get(key).getDistances().get(indexOfGeneThatNeedsDistances).get(indexOfDistanceValue));
//                }
//            }
//        }
//
//        allSpeciesOut.sort(Comparator.comparing(Species::getName));
//        allSpeciesOut.forEach(s -> s.getGenes().sort(Comparator.comparing(Gene::getId)));
//        log.warning("Running for " + outFileName);
//        List<GeneTuple> foundHgtsOut = DistributionAlgorithm.runAlgorithm(allSpeciesOut, percentage, multithreading, runMode);
//        NexusWriter.writeHgtsToFile("tmp", nrOut.getDistanceMap(), foundHgtsOut);
//
//        if (outputfilename != null) {
//            NexusWriter.writeHgtsToFile(outputfilename, nr.getDistanceMap(), foundHgts);
//            log.setLevel(Level.WARNING);
//            Utils.evaluateResults(allSpeciesOut, nrOut.getGenesWithRCBRelations());
//
//            Utils.evaluateResults(allSpecies, nrOut.getGenesWithRCBRelations());
//        }

        log.warning("Tool finished.");

        log.warning("Time elapsed: " + ((System.currentTimeMillis() - start) / 1000F) + " seconds.");
    }
}
