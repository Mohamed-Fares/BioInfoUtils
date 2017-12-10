package edu.nu.bio.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;

import org.apache.commons.lang3.StringUtils;

public class FastaToFastqKmers {

	public static void main(String[] args) {

		generateSimulatedReads(args);

	}

	private static void generateSimulatedReads(String[] args) {
		String pathName;
		if(args.length < 1){
			System.out.println("too few arguments !!!");
//			System.out.println("Enter the input fasta file name ");
			
//			Debug
//			String FileName = "Homo_sapiens.GRCh38.cds.all.fa";
//			pathName = "D:\\Research\\Nile_University\\Research\\RNA_seq_assembly\\Khmer_assessment\\" + FileName;
			
//			String FileName = "Homo_sapiens.GRCh38.cds.all.fa";
			String FileName = "HS_GRCh38.cdna_ncrna_comb.fa";
			pathName = "D:\\Research\\Nile_University\\Research\\RNA_seq_assembly\\Khmer_assessment\\fastq_seq_concat\\" + FileName;
		} else {
			pathName = args[0];
			
			
		}
		
		File file = new File(pathName);
		System.out.println("input file : " + file.getAbsolutePath());

		BufferedReader reader = null;
		PrintWriter writer = null;
		PrintWriter logger = null;

		try {
			reader = new BufferedReader(new FileReader(file));
			writer = new PrintWriter(new FileWriter(pathName + ".fastq.test"));
			logger = new PrintWriter(new FileWriter(pathName + ".log"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		String transcriptName;
		String genetName;
		int transcriptNumber = 0;
		long readsCount = 0;
		HashSet geneNamesSet = new HashSet();
		
		String line;

		try {
			line = reader.readLine();
			do {
				if (line.startsWith(">")) {
					transcriptName = StringUtils.substringBefore(line, " ")
							.substring(1);
					
					genetName = StringUtils.substringBetween(line, "gene:", " ");
					
					geneNamesSet.add(transcriptName);
					
					line = reader.readLine();
					StringBuilder seqBuilder = new StringBuilder();
					do {
						seqBuilder.append(line);
						line = reader.readLine();
					} while (line != null && !line.startsWith(">"));

					String[] reads = generateReadsFromTranscript(
							seqBuilder.toString(), 100);
					for (int i = 0; i < reads.length; i++) {
//						printReadSeq(writer, genetName + "_" + transcriptName + "_"+ transcriptNumber + "_" + i, reads[i]);
						printFasaReadSeq(writer, "g:" + genetName + "_t:" + transcriptName + "_i:" +i, reads[i]);
					}
					
					System.out.println("generating transcript : "+ transcriptNumber);

					readsCount += reads.length;
					logger.print("gene : " + transcriptName);
					logger.print(" , gene number : " + transcriptNumber);
					logger.print(" , Number of reads : "+ reads.length);
					logger.print(" , size of transcript : "+ seqBuilder.length());
					logger.println();
					logger.println("Total Reads count : " + readsCount);

					transcriptNumber++;
					writer.flush();
					logger.flush();
				} else {
					line = reader.readLine();
				}
			} while (line != null);
			
			
			logger.println("Total genes count : " + geneNamesSet.size());
			logger.println("Total Reads count : " + readsCount);
			writer.flush();
			logger.flush();

			reader.close();
			writer.close();
			logger.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String[] generateReadsFromTranscript(String trans,
			int readlength) {
		String[] reads = null;
		if (trans.length() <= readlength) {
//			reads = new String[] { trans };
			reads = new String[0];
		} else {
			int readsNumber = trans.length() - readlength + 1;
			reads = new String[readsNumber];
			for (int i = 0; i < readsNumber; i++) {
				reads[i] = StringUtils.substring(trans, i, i + readlength);
			}
		}
		return reads;
	}

	private static void printReadSeq(PrintWriter writer, String id, String seq) {
		writer.print("@");
		writer.print(id);
		writer.println();
		writer.println(seq);
		writer.println("+");
		writer.println(createFakeQualityScores(seq));
	}

	private static void printFasaReadSeq(PrintWriter writer, String id, String seq) {
		writer.print(">");
		writer.print(id);
		writer.println();
		writer.print(seq);
		writer.println();
		
	}

	private static String createFakeQualityScores(String seq) {
		return seq.replaceAll(".", "H");
	}



}
