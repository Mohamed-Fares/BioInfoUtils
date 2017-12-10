package edu.nu.bio.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

public class IterateKhmerSpanKmerSizes {

	public static void main(String[] args) {

		if (args.length > 0) {
			final String inputFile = args[0];
			new Thread() {
				public void run() {
					for (int i = 17; i < 76; i+=2) {
						runKhmer(i, inputFile);
						summarizeResults(inputFile+".part", inputFile+"_summary_"+i);
					}

					
//					for (int i = 17; i < 18; i += 2) {
//						runKhmer(i, inputFile);
//						summarizeResults(inputFile+".part", inputFile+"_summary_"+i);
//						
//					}
				}
			}.start();
		} else {
			System.err.println("Too few arguments ..... ");
		}
		
//		String path = "D:\\Research\\Nile_University\\Research\\RNA_seq_assembly\\Khmer_assessment\\test\\";
//		summarizeResults(path+"test.fa.part", path+"output");

	}

	private static String runKhmer(int kmerSize, String inputFile) {

		Runtime rt = Runtime.getRuntime();
		String command1 = "do-partition.py -k " + kmerSize
				+ " -T 4 -M 30000000000 --force output_"+ kmerSize+ " " + inputFile;
		System.out.println(command1);
		Process pr;
		try {
			pr = rt.exec(command1);
			
			BufferedReader stdInput = new BufferedReader(new InputStreamReader(pr.getInputStream()));

			BufferedReader stdError = new BufferedReader(new InputStreamReader(pr.getErrorStream()));

				// read the output from the command
			
			try {
				pr.waitFor();
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
			
			
			PrintWriter log = new PrintWriter(inputFile+"_"+kmerSize+".log");
//			System.out.println("Here is the standard output of the command:\n");
			log.write("Here is the standard output of the command:\n");
			String s = null;
			while ((s = stdInput.readLine()) != null) {
				System.out.println(s+"\n");
				log.write(s+"\n");
			}
			log.flush();

			
			log.write("Here is the standard error of the command (if any):\n");
			while ((s = stdError.readLine()) != null) {
				System.err.println(s+"\n");
				log.write(s+"\n");
			}
			log.flush();
			log.close();
			
			while (pr.isAlive()) {
				try {
					Thread.currentThread().sleep(5000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		return inputFile + ".part";

	}

	private static void summarizeResults(String inputFile, String outputFile) {


		try {
			File file = new File(inputFile);
			int counter = 0;
			while (!file.exists() && counter < 20) {
				System.err.println("file doesn't exists .... ");
				System.err.println("waiting for thread to finish.... ");
				try {
					Thread.currentThread().sleep(5000);
					
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				counter++;
			}
			BufferedReader reader = new BufferedReader(new FileReader(file));
			PrintWriter writer = new PrintWriter(new File(outputFile));
			System.out.println("reading partitioned file ..........");

			HashSet<TranscriptInfo> nodesSet = new HashSet<>();
			String gene_name;
			String transcript_id;
			int gene_cluster;

			String line;
			int progress = 0;
			try {
				line = reader.readLine();
				do {
					if (line.startsWith(">")) {

						gene_name = StringUtils.substringBetween(line, "g:",
								"_t:");

						transcript_id = StringUtils.substringBetween(line,
								"_t:", "_i:");

						String cluster = StringUtils.substringAfter(line, "	")
								.trim();
						gene_cluster = Integer.parseInt(cluster);

						TranscriptInfo info = new TranscriptInfo();
						info.geneID = gene_name;
						info.transcriptId = transcript_id;
						info.cluster = gene_cluster;
						nodesSet.add(info);

						reader.readLine();
						line = reader.readLine();

						progress++;
						if (progress % 100 == 0) {
							System.out.println("progress : " + progress);
						}

					} else {
						line = reader.readLine();
					}
				} while (line != null);

				reader.close();
				file.delete();
				System.out.println("reading source finsished ..........");

				progress = 0;

				System.out.println("writing summary ......");

				for (TranscriptInfo transcriptInfo : nodesSet) {
					writer.print(transcriptInfo.transcriptId);
					writer.print("," + transcriptInfo.geneID);
					writer.print("," + transcriptInfo.cluster);
					writer.println();
				}

				writer.flush();

				writer.close();
				
				
				
				generateStats(nodesSet, outputFile);

			} catch (IOException e) {
				e.printStackTrace();
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}
	
	
	private static void generateStats(HashSet<TranscriptInfo> nodes, String fileName){
		
//		File parentDirectory = new File(fileName).getParentFile();
//		
//		File genesPerPartitionDir = new File(parentDirectory.getAbsolutePath()+File.pathSeparator+"gpp");
//		if(!genesPerPartitionDir.exists()){
//			genesPerPartitionDir.mkdir();
//		}
//		File transcriptsPerPartitionDir = new File(parentDirectory.getAbsolutePath()+File.pathSeparator+"tpp");
//		if(!transcriptsPerPartitionDir.exists()){
//			transcriptsPerPartitionDir.mkdir();
//		}
//		File partitionsPerGeneDir = new File(parentDirectory.getAbsolutePath()+File.pathSeparator+"ppg");
//		if(!partitionsPerGeneDir.exists()){
//			partitionsPerGeneDir.mkdir();
//		}
//		File partitionsPerTranscriptDir = new File(parentDirectory.getAbsolutePath()+File.pathSeparator+"ppt");
//		if(!partitionsPerTranscriptDir.exists()){
//			partitionsPerTranscriptDir.mkdir();
//		}
		
		Map<Integer, HashSet<String>> genesPerPartition = new HashMap<Integer, HashSet<String>>();
		Map<Integer, HashSet<String>> transcriptsPerPartition = new HashMap<Integer, HashSet<String>>();
		Map<String, HashSet<Integer>> partitionsPerGene = new HashMap<String, HashSet<Integer>>();
		Map<String, HashSet<Integer>> partitionsPerTranscript = new HashMap<String, HashSet<Integer>>();
		
		
		for (TranscriptInfo info : nodes) {
			add(transcriptsPerPartition, info.cluster, info.transcriptId);
			add(genesPerPartition, info.cluster, info.geneID);
			add(partitionsPerTranscript, info.transcriptId, info.cluster);
			add(partitionsPerGene, info.geneID, info.cluster);
		}
		
		
		
		try {
			PrintWriter genesPerPartitionWriter = new PrintWriter(new File(fileName+"_gpp"+".csv"));
			PrintWriter transcriptsPerPartitionWriter = new PrintWriter(new File(fileName+"_tpp"+".csv"));
			PrintWriter partitionPerGeneWriter = new PrintWriter(new File(fileName+"_ppg"+".csv"));
			PrintWriter partitionPerTranscriptWriter = new PrintWriter(new File(fileName+"_ppt"+".csv"));

//			PrintWriter genesPerPartitionWriter = new PrintWriter(new File(genesPerPartitionDir.getAbsoluteFile()+File.separator+fileName+".csv"));
//			PrintWriter transcriptsPerPartitionWriter = new PrintWriter(new File(transcriptsPerPartitionDir.getAbsoluteFile()+File.separator+fileName+".csv"));
//			PrintWriter partitionPerGeneWriter = new PrintWriter(new File(partitionsPerGeneDir.getAbsoluteFile()+File.separator+fileName+".csv"));
//			PrintWriter partitionPerTranscriptWriter = new PrintWriter(new File(partitionsPerTranscriptDir.getAbsoluteFile()+File.separator+fileName+".csv"));
			
			
			
			
			writeMapToFile((HashMap<Integer, HashSet<String>>)genesPerPartition, genesPerPartitionWriter);
			writeMapToFile((HashMap<Integer, HashSet<String>>)transcriptsPerPartition, transcriptsPerPartitionWriter);
			
			
			
			writeMapToFile(partitionPerGeneWriter, (HashMap<String, HashSet<Integer>>)partitionsPerGene);
			writeMapToFile(partitionPerTranscriptWriter, (HashMap<String, HashSet<Integer>>)partitionsPerTranscript);
			
			genesPerPartitionWriter.flush();
			transcriptsPerPartitionWriter.flush();
			partitionPerGeneWriter.flush();
			partitionPerTranscriptWriter.flush();

			genesPerPartitionWriter.close();
			transcriptsPerPartitionWriter.close();
			partitionPerGeneWriter.close();
			partitionPerTranscriptWriter.close();
			
			
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		
		
	}
	
	
	private static void writeMapToFile(PrintWriter writer, HashMap<String, HashSet<Integer>> map) {
		Set<String> keys = map.keySet();
		for (String key : keys) {
			HashSet<Integer> values = (HashSet<Integer>) map.get(key);
			writer.print(key);
			writer.print("," + values.size());
			writer.println();
		}

	}

	private static void writeMapToFile(HashMap<Integer, HashSet<String>> map, PrintWriter writer) {
		Set<Integer> keys = map.keySet();
		for (Integer key : keys) {
			HashSet<String> values = (HashSet<String>) map.get(key);
			writer.print(key);
			writer.print("," + values.size());
			writer.println();
		}
		
	}
	
	private static void add(Map<String, HashSet<Integer>> container, String key, Integer value) {
	    HashSet<Integer> values = container.get(key);
	    if (values == null) {
	        values = new HashSet<Integer>();
	    }
	    values.add(value);
	    container.put(key, values);
	}
	
	
	private static void add(Map<Integer, HashSet<String>> container, Integer key, String value) {
	    HashSet<String> values = container.get(key);
	    if (values == null) {
	        values = new HashSet<String>();
	    }
	    values.add(value);
	    container.put(key, values);
	}

}
