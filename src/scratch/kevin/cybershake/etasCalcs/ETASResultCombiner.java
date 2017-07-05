package scratch.kevin.cybershake.etasCalcs;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.Collections;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.util.FileUtils;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;

import com.google.common.base.Preconditions;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.io.CharStreams;
import com.google.common.io.Files;

public class ETASResultCombiner {

	public static void main(String[] args) throws IOException {
		// used to combine multiple etas result zip files into a single zip file
		
		File simsDir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/sims");
//		String prefix = "2015_06_15-u2mapped-bombay_beach_brawley_fault_m6";
//		String prefix = "2015_06_15-u2mapped-parkfield";
		String prefix = "2015_06_15-u2mapped-mojave_s_point_m6";
		int[] rounds = { 1, 2 };
		
		File[] zipFiles = new File[rounds.length];
		
		for (int i=0; i<rounds.length; i++)
			zipFiles[i] = new File(new File(simsDir, prefix+"-round"+rounds[i]), "results.zip");
		
		
		File outputFile = new File(simsDir, prefix+"-combined.zip");
		
		boolean onlyWithFaultBased = false;
		double minMag = 5d;
		
//		boolean onlyWithFaultBased = true;
//		double minMag = 0d;
		
		if (onlyWithFaultBased) {
			String name = outputFile.getName();
			Preconditions.checkState(name.endsWith(".zip"));
			name = name.substring(0, name.indexOf(".zip"))+"_only_with_fault.zip";
			outputFile = new File(outputFile.getParentFile(), name);
		}
		
		File workDir = Files.createTempDir();
		
		List<List<ETAS_EqkRupture>> catalogs = Lists.newArrayList();
		
		List<String> infoStrings = Lists.newArrayList();
		
		for (File zipFile : zipFiles) {
			System.out.print("Loading "+zipFile.getAbsolutePath()+" ... ");
			ZipFile zip = new ZipFile(zipFile);
			
			for (ZipEntry entry : Collections.list(zip.entries())) {
				String name = entry.getName();
//				if (Math.random() < 0.0001)
//					System.out.println(name);
				if (!entry.isDirectory()) {
					// fix for ones where the slash was left out accidentally on directory entires
					if (name.startsWith("sim_") && Character.isDigit(name.charAt(name.length()-1)))
						name = name+"/";
					else
						continue;
				}
//				if (Math.random() > 0.01)
//					continue;
//				System.out.println(entry.getName());
				String subEntryName = name+"simulatedEvents.txt";
				ZipEntry catEntry = zip.getEntry(subEntryName);
				String infoEntryName = name+"infoString.txt";
				ZipEntry infoEntry = zip.getEntry(infoEntryName);
				if (catEntry == null || infoEntry == null)
					continue;
				
				// make sure it's actually done
				BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(infoEntry)));
				
				StringBuilder infoString = new StringBuilder();
				boolean done = false;
				for (String line : CharStreams.readLines(reader)) {
					infoString.append(line);
					infoString.append("\n");
					if (line.contains("Total num ruptures: ")) {
						done = true;
						break;
					}
				}
				if (!done)
					continue;
//				System.out.println("Loading "+catEntry.getName());
				
				List<ETAS_EqkRupture> catalog;
				try {
					catalog = ETAS_CatalogIO.loadCatalog(zip.getInputStream(catEntry), minMag);
				} catch (Exception e) {
					continue;
				}
				// now remove all spontaneous
				catalog = ETAS_SimAnalysisTools.getChildrenFromCatalog(catalog, 0);
				
				if (onlyWithFaultBased) {
					boolean found = false;
					for (ETAS_EqkRupture rup : catalog) {
						if (rup.getFSSIndex() >= 0) {
							found = true;
							break;
						}
					}
					if (!found)
						continue;
				}
				
				catalogs.add(catalog);
				infoStrings.add(infoString.toString());
			}
			System.out.println(catalogs.size()+" catalogs.");
			zip.close();
		}
		
		System.out.println("Loaded "+catalogs.size()+" catalogs");
		
		// now write new ones
		int numDigits = (""+(catalogs.size()-1)).length();
		
		List<String> zipNames = Lists.newArrayList();
		
		for (int i=0; i<catalogs.size(); i++) {
			String name = i+"";
			while (name.length() < numDigits)
				name = "0"+name;
			name = "sim_"+name;
			
			File resultsDir = new File(workDir, name);
			resultsDir.mkdir();
			File catFile = new File(resultsDir, "simulatedEvents.txt");
			File infoFile = new File(resultsDir, "infoString.txt");
			ETAS_CatalogIO.writeEventDataToFile(catFile, catalogs.get(i));
			Files.write(infoStrings.get(i), infoFile, Charset.defaultCharset());
			zipNames.add("/"+name+"/");
			zipNames.add("/"+name+"/simulatedEvents.txt");
			zipNames.add("/"+name+"/infoString.txt");
			resultsDir.mkdir();
		}
		
		FileUtils.createZipFile(outputFile.getAbsolutePath(), workDir.getAbsolutePath()+"/", zipNames);
		
		FileUtils.deleteRecursive(workDir);
	}

}
