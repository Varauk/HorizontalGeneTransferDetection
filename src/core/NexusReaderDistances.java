package core;

import objects.Gene;
import objects.GeneDistances;
import objects.GeneTuple;
import objects.Matrix;
import objects.NexusFormatException;
import objects.Species;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NexusReaderDistances {

    private HashMap<String, GeneDistances> distanceMap;

    private List<Gene> genesWithRCBRelations = new ArrayList<>();

    /**
     * Opens a nexus file and reads in all data that is important to this program.
     *
     * @param filePath path to the nexus file
     * @throws IOException          error reading the file
     * @throws NexusFormatException error parsing the file
     */
    public NexusReaderDistances(String filePath) throws IOException, NexusFormatException {

        BufferedReader bufferedReader = null;
        FileReader fileReader = null;
        try {
            fileReader = new FileReader(filePath);
            bufferedReader = new BufferedReader(fileReader);

            int state = 0; //0=start of file, 1=after #NEXUS. 2=inside block
            boolean commentMode = false;
            TreeMap<Integer, String> blockLines = null;
            String blockName = null;
            Matcher beginBlockMatcher = Pattern.compile("^(?i)BEGIN ([A-Z]+);").matcher("");
            String line = bufferedReader.readLine();
            int lineNumber = 0;
            while (line != null) {
                lineNumber++;
                line = line.trim();
                //end multiline comments
                if (commentMode) {
                    int endIndex = line.indexOf("]");
                    if (endIndex == -1) {
                        line = bufferedReader.readLine();
                        continue;
                    } else {
                        line = line.substring(endIndex);
                        commentMode = false;
                    }
                }
                //remove one line comments
                line = line.replaceAll("\\[.*]", "");
                //remove multi line comments
                int startIndex = line.indexOf("[");
                if (startIndex != -1) {
                    line = line.substring(0, startIndex);
                    commentMode = true;
                }
                //skip emtpy lines
                if (line.isEmpty()) {
                    line = bufferedReader.readLine();
                    continue;
                }

                //===read lines===
                switch (state) { //0=start of file, 1=after #NEXUS. 2=inside block
                    case 0:
                        if (!line.equalsIgnoreCase("#NEXUS")) {
                            throw new NexusFormatException("File must start with #NEXUS. Line: " + lineNumber);
                        } else {
                            state++;
                        }
                        break;
                    case 1:
                        beginBlockMatcher.reset(line);
                        if (beginBlockMatcher.matches()) {
                            blockLines = new TreeMap<Integer, String>();
                            blockName = beginBlockMatcher.group(1);
                            state++;
                        } else {
                            throw new NexusFormatException("Error parsing BEGIN statement. Line " + lineNumber);
                        }
                        break;
                    case 2:
                        if (line.equalsIgnoreCase("END;")) {
                            makeBlock(blockName, blockLines);
                            state--;
                        } else {
                            blockLines.put(lineNumber, line);
                        }
                        break;
                    default:
                        throw new RuntimeException("NexusReader reached illegal state: " + state);
                }

                line = bufferedReader.readLine();
            }
            if (state == 2) {
                throw new NexusFormatException("Didn't close block at end of file.");
            }
            if (state == 0) {
                throw new NexusFormatException("File seems empty.");
            }

            bufferedReader.close();

        } catch (IOException e) {
            if (bufferedReader != null) {
                bufferedReader.close();
            }
            if (fileReader != null) {
                fileReader.close();
            }
            throw e;
        } catch (NexusFormatException e) {
            if (bufferedReader != null) {
                bufferedReader.close();
            }
            if (fileReader != null) {
                fileReader.close();
            }
            throw e;
        }
    }


    private void makeBlock(String blockName, Map<Integer, String> lines)
            throws NexusFormatException {
        if (blockName.equalsIgnoreCase("ALLDISTANCES")) {
            readInAllDistances(lines);
        } else if (blockName.equalsIgnoreCase("RELATIONS")) {
            readInRelations(lines);
        }
        //ignore unknown blocks
    }

    private void readInAllDistances(Map<Integer, String> lines) throws NexusFormatException {
        Logger log = Logger.getGlobal();
        int state = 0; //0=out of command, 1=in config, 2=in matrix, 3=end of matrix
        if (distanceMap == null) {
            distanceMap = new HashMap<String, GeneDistances>();
        }

        String curName = null;
        Matrix.MatrixFormat curFormat = null;
        GeneDistances curDistances = null;

        for (Map.Entry<Integer, String> lineEntry : lines.entrySet()) {
            int lineNumber = lineEntry.getKey();
            String line = lineEntry.getValue();

            String[] commands = line.split("\\s+");
            if (state == 0 && commands[0].equalsIgnoreCase("distances")) {
                state++;
                if (commands.length > 1) {
                    line = line.substring(9);
                    line = line.trim();
                    commands = line.split("\\s+");
                } else {
                    continue;
                }
            }
            if (state == 1) {
                for (String command : commands) {
                    String[] data = command.split("=");
                    if (data.length != 2) {
                        throw new NexusFormatException("No data given to sub-command. Line: " + lineNumber);
                    }
                    if (data[0].equalsIgnoreCase("name")) {
                        curName = data[1];
                    } else if (data[0].equalsIgnoreCase("triangle")) {
                        try {
                            curFormat = Matrix.MatrixFormat.valueOf(data[1].toUpperCase());
                        } catch (IllegalArgumentException e) {
                            throw new NexusFormatException(
                                    String.format("Matrix format unknown: %s. Line: %d", data[1], lineNumber));
                        }
                    }
                    //ignore unknown sub-commands
                }

                if (curName != null && curFormat != null) {
                    curDistances = new GeneDistances(curFormat, curName);
                    state++;
                }
            } else if (state == 2) {
                //read in matrix
                if (line.equals(";")) {
                    state++;
                } else {
                    String[] speciesAndGeneName = commands[0].split("/");
                    Species species = curDistances.getSpecies().stream().filter(speciesL -> speciesL.getName().equals(speciesAndGeneName[0])).findFirst().orElse(null);
                    if (species == null) {
                        curDistances.getSpecies().add(new Species(speciesAndGeneName[0]));
                        species = curDistances.getSpecies().stream().filter(speciesL -> speciesL.getName().equals(speciesAndGeneName[0])).findFirst().orElse(null);
                    }
                    if (speciesAndGeneName.length == 1) {
                        log.severe("Malformed entry in input file! Gene tree: " + curName + ", line number: " + lineNumber);
                        System.exit(1);
                    }
                    curDistances.getGenes().add(new Gene(speciesAndGeneName[1], species, curName));
                    ArrayList<Double> row = new ArrayList<Double>();
                    for (int i = 1; i < commands.length; i++) {
                        if (commands[i].equals(";")) {
                            continue;
                        } else if (commands[i].endsWith(";")) {
                            commands[i] = commands[i].replace(";", "");
                        }
                        try {
                            row.add(Double.valueOf(commands[i]));
                        } catch (NumberFormatException e) {
                            throw new NexusFormatException("Cannot read double number. Line: " + lineNumber, e);
                        }
                    }
                    curDistances.getDistances().add(row);

                    if (line.endsWith(";")) {
                        state++;
                    }
                }
                if (state == 3) {
                    //add distances to map
                    distanceMap.put(curName, curDistances);
                    //accept new ones
                    state = 0;
                }
            }
            //ignore other commands
        }
    }

    private void readInRelations(Map<Integer, String> lines)
            throws NexusFormatException {

        int state = 0; //0=out of command, 1=in config, 2=in matrix, 3=end of matrix
        //relationMap = new HashMap<GeneRelationType, List<GeneRelation>>();
        boolean rcb = false;

        //GeneRelationType curType = null;
        String curName = null;
        Matrix.MatrixFormat curFormat = null;

        List<Gene> genesInCurrentTree = null;

        //GeneRelation curRelation = null;
        for (Map.Entry<Integer, String> lineEntry : lines.entrySet()) {

            int lineNumber = lineEntry.getKey();
            String line = lineEntry.getValue();

            String[] commands = line.split("\\s+");
            if (state == 0 && commands[0].equalsIgnoreCase("relation")) {
                state++;
                if (commands.length > 1) {
                    line = line.substring(9);
                    line = line.trim();
                    commands = line.split("\\s+");
                } else {
                    continue;
                }
            }
            if (state == 1) {

                for (String command : commands) {
                    String[] data = command.split("=");
                    if (data.length != 2) {
                        throw new NexusFormatException("No data given to sub-command. Line: " + lineNumber);
                    }
                    if (data[0].equalsIgnoreCase("type")) {
                        if (data[1].equalsIgnoreCase("rcb")) {
                            rcb = true;
                        } else {
                            rcb = false;
                        }

                    } else if (data[0].equalsIgnoreCase("name")) {
                        curName = data[1];
                        if (genesInCurrentTree != null) {
                            for (int i = 0; i < genesInCurrentTree.size(); i++) {
                                for (int j = 1; j < genesInCurrentTree.get(i).getRelationsToOtherGenesInItsTreeAsString().length; j++) {
                                    GeneTuple.RelationType type = GeneTuple.RelationType.EQUAL;
                                    if (genesInCurrentTree.get(i).getRelationsToOtherGenesInItsTreeAsString()[j].equals("-1")) {
                                        type = GeneTuple.RelationType.LESSER;
                                    } else if (genesInCurrentTree.get(i).getRelationsToOtherGenesInItsTreeAsString()[j].equals("1")) {
                                        type = GeneTuple.RelationType.HIGHER;
                                    }
                                    genesInCurrentTree.get(i).getRelationsToOtherGenesInItsTree().add(new GeneTuple(genesInCurrentTree.get(i), genesInCurrentTree.get(j - 1), type));
                                }
                            }
                        }
                        genesInCurrentTree = new ArrayList<>();
                    } else if (data[0].equalsIgnoreCase("triangle")) {
                        try {
                            curFormat = Matrix.MatrixFormat.valueOf(data[1].toUpperCase());
                        } catch (IllegalArgumentException e) {
                            throw new NexusFormatException(
                                    String.format("Matrix format unknown: %s. Line: %d", data[1], lineNumber));
                        }
                    }
                    //ignore unknown sub-commands
                }

                if (curName != null && curFormat != null) {
                    // = new GeneRelation(curType, curFormat, curName);
                    state++;
                }
            } else if (state == 2) {
                //read in matrix
                if (line.equals(";")) {
                    state++;
                } else {
                    Species tmp = new Species(commands[0].split("/")[0]);
                    Gene gene = new Gene(commands[0].split("/")[1], tmp, curName);
                    for (int i = 1; i < commands.length; i++) {
                        if (commands[i].equals(";")) {
                            continue;
                        } else if (commands[i].endsWith(";")) {
                            commands[i] = commands[i].replace(";", "");
                        }
                        gene.setRelationsToOtherGenesInItsTreeAsString(commands);
                    }
                    if (rcb) {
                        genesWithRCBRelations.add(gene);
                        genesInCurrentTree.add(gene);
                    }

                    //curRelation.getRelation().add(byteList);

                    if (line.endsWith(";")) {
                        state++;
                    }
                }
                if (state == 3) {
                    //add relation to list
                    /*List<GeneRelation> relationList = relationMap.get(curType);
                    if (relationList == null) {
                        relationList = new ArrayList<GeneRelation>();
                        relationMap.put(curType, relationList);
                    }
                    relationList.add(curRelation);*/
                    //accept new ones
                    state = 0;
                }
            }
            //ignore other commands
        }

    }

    public HashMap<String, GeneDistances> getDistanceMap() {
        return distanceMap;
    }

    public List<Gene> getGenesWithRCBRelations() {
        return genesWithRCBRelations;
    }
}
