package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationPageGen.ValidationResult;

public class BBP_PartBSummaryPageGen {
	
	private static class ResultKey {
		private final VelocityModel vm;
		private final boolean isRotation;
		private final Scenario scenario;
		private final float distance;
		
		public ResultKey(VelocityModel vm, boolean isRotation, Scenario scenario, float distance) {
			super();
			this.vm = vm;
			this.isRotation = isRotation;
			this.scenario = scenario;
			this.distance = distance;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Float.floatToIntBits(distance);
			result = prime * result + (isRotation ? 1231 : 1237);
			result = prime * result + ((scenario == null) ? 0 : scenario.hashCode());
			result = prime * result + ((vm == null) ? 0 : vm.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			ResultKey other = (ResultKey) obj;
			if (Float.floatToIntBits(distance) != Float.floatToIntBits(other.distance))
				return false;
			if (isRotation != other.isRotation)
				return false;
			if (scenario != other.scenario)
				return false;
			if (vm != other.vm)
				return false;
			return true;
		}
	}
	
	private static String bbpTableStr(boolean passes, boolean official) {
		if (official) {
			if (passes)
				return "**PASS**";
			return "**FAIL**";
		}
		if (passes)
			return "*(PASS)*";
		return "*(FAIL)*";
	}
	
	private static boolean isOfficial(Scenario scenario, float distance) {
		return scenario.isOfficialCriteria() && isOfficialPartB_Distance(distance);
	}
	
	private static boolean isOfficialPartB_Distance(float distance) {
		for (double d : BBP_PartBValidationConfig.OFFICIAL_DISTANCES)
			if ((float)d == distance)
				return true;
		return false;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_3165.instance(baseDir);
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		Map<ResultKey, File> resultDirMap = new HashMap<>();
		Map<ResultKey, List<ValidationResult>> resultsMap = new HashMap<>();
		Map<ResultKey, String> plotsMap = new HashMap<>();
		
		// look for result
		for (VelocityModel vm : VelocityModel.values()) {
			File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
			if (!vmDir.exists())
				continue;
			System.out.println("Processing "+vm+" from "+vmDir.getAbsolutePath());
			File partBdir = new File(vmDir, "bbp_part_b");
			if (partBdir.exists()) {
				File resourcesDir = new File(partBdir, "resources");
				File csvFile = new File(resourcesDir, "part_b_results.csv");
				if (csvFile.exists()) {
					System.out.println("Found regular PartB CSV file: "+csvFile.getAbsolutePath());
					Table<Scenario, Double, String> plotTable = HashBasedTable.create();
					Table<Scenario, Double, List<ValidationResult>> resultTable = BBP_PartBValidationPageGen.loadResultsCSV(csvFile, plotTable);
					for (Cell<Scenario, Double, List<ValidationResult>> cell : resultTable.cellSet()) {
						Scenario scenario = cell.getRowKey();
						Double distance = cell.getColumnKey();
						List<ValidationResult> results = cell.getValue();
						String plotName = plotTable.get(scenario, distance);
						ResultKey key = new ResultKey(vm, false, scenario, distance.floatValue());
						resultDirMap.put(key, partBdir);
						resultsMap.put(key, results);
						plotsMap.put(key, plotName);
					}
				}
			}
			// now look for rotation results
			for (Scenario scenario : Scenario.values()) {
				File rotationDir = new File(vmDir, "rotated_ruptures_"+scenario.getPrefix());
				if (rotationDir.exists()) {
					File resourcesDir = new File(rotationDir, "resources");
					File csvFile = new File(resourcesDir, "part_b_results.csv");
					if (csvFile.exists()) {
						System.out.println("Found rotation PartB CSV file: "+csvFile.getAbsolutePath());
						Table<Scenario, Double, String> plotTable = HashBasedTable.create();
						Table<Scenario, Double, List<ValidationResult>> resultTable = BBP_PartBValidationPageGen.loadResultsCSV(csvFile, plotTable);
						for (Cell<Scenario, Double, List<ValidationResult>> cell : resultTable.cellSet()) {
							Preconditions.checkState(scenario == cell.getRowKey());
							Double distance = cell.getColumnKey();
							List<ValidationResult> results = cell.getValue();
							String plotName = plotTable.get(scenario, distance);
							ResultKey key = new ResultKey(vm, true, scenario, distance.floatValue());
							resultDirMap.put(key, rotationDir);
							resultsMap.put(key, results);
							plotsMap.put(key, plotName);
						}
					}
				}
			}
		}
		
		System.out.println("Found "+resultDirMap.size()+" studies");
		
		if (resultDirMap.isEmpty())
			System.exit(2);
		
		File summaryDir = new File(catalogOutputDir, "bbp_part_b_summary");
		Preconditions.checkState(summaryDir.exists() || summaryDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# "+catalog.getName()+" BBP PartB Results Summary");
		lines.add("");
		lines.add("This page combines results from multiple BBP PartB calculations for different BBP velocity models and rupture "
				+ "rotation studies. Refer to individual pages for details and methodology for individual calculations. Regular "
				+ "PartB studies distribute sites around RSQSim ruptures in a half-racetrack. Rotation studies instead use a fixed "
				+ "site location, translating each rupture to the specified distance from the site, and rotating the ruptures about "
				+ "their moment centroids to sample many source azimuths.");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		lines.add("## BBP PartB Background Information");
		lines.add(topLink); lines.add("");
		lines.addAll(BBP_PartBValidationPageGen.getBackgroundLines());
		lines.add("");
		
		HashSet<Float> distancesSet = new HashSet<>();
		HashSet<Scenario> scenariosSet = new HashSet<>();
		HashSet<VelocityModel> vmsSet = new HashSet<>();
		for (ResultKey key : resultsMap.keySet()) {
			distancesSet.add(key.distance);
			scenariosSet.add(key.scenario);
			vmsSet.add(key.vm);
		}
		List<Float> distances = new ArrayList<>(distancesSet);
		Collections.sort(distances);
		
		lines.add("## Result");
		lines.add(topLink); lines.add("");
		lines.add("Results for official BBP PartB criteria are listed in **bold**, and those for unofficial scenarios (which use the same "
				+ "formulae to determine criteria but where the underlying models are less constrained) are listed in *(italics)*. Failures "
				+ "also list the largetst failure (at any period), in natural-log units from the criterion. A positive value means that the "
				+ "simulated median value was above the maximum criterion by the specified natural-log amount, and a negative below the minimum "
				+ "criterion.");
		lines.add("");
		
		for (Scenario scenario : Scenario.values()) {
			if (!scenariosSet.contains(scenario))
				continue;
			
			lines.add("### "+scenario.getName());
			lines.add(topLink); lines.add("");
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("Velocity Model");
			table.addColumn("Rotation?");
			table.addColumn("Link");
			for (Float distance : distances)
				table.addColumn(distance+" km");
			table.finalizeLine();
			
			for (VelocityModel vm : VelocityModel.values()) {
				if (!vmsSet.contains(vm))
					continue;
				for (boolean rotation : new boolean[] { false, true }) {
					File resultDir = null;
					for (Float distance : distances) {
						ResultKey key = new ResultKey(vm, rotation, scenario, distance);
						if (resultsMap.containsKey(key)) {
							resultDir = resultDirMap.get(key);
							break;
						}
					}
					if (resultDir != null) {
						table.initNewLine();
						table.addColumn(vm);
						table.addColumn(rotation ? "yes" : "no");
						String link = "../"+resultDir.getParentFile().getName()+"/"+resultDir.getName();
						table.addColumn("[Page Link]("+link+")");
						List<String> imageLinks = new ArrayList<>();
						for (Float distance : distances) {
							ResultKey key = new ResultKey(vm, rotation, scenario, distance);
							if (resultsMap.containsKey(key)) {
								boolean passes = true;
								double maxLnDiff = 0;
								double maxPeriod = Double.NaN;
								for (ValidationResult result : resultsMap.get(key)) {
									passes = passes && result.passes();
									if (!result.passes()) {
										double lnDiff = result.getLogFailAmount();
										if (Math.abs(lnDiff) > Math.abs(maxLnDiff)) {
											maxLnDiff = lnDiff;
											maxPeriod = result.getPeriod();
										}
									}
								}
								String passStr = bbpTableStr(passes, isOfficial(scenario, distance));
								if (!passes)
									passStr +=", Max Ln Fail: "+(float)maxLnDiff+" @ "+(float)maxPeriod+"s";
								table.addColumn(passStr);
								imageLinks.add("![Plot]("+link+"/resources/"+plotsMap.get(key)+")");
							} else {
								table.addColumn("*N/A*");
								imageLinks.add("");
							}
						}
						table.finalizeLine();
						table.initNewLine();
						table.addColumn("");
						table.addColumn("");
						table.addColumn("");
						for (String image : imageLinks)
							table.addColumn(image);
						table.finalizeLine();
					}
				}
			}
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, summaryDir);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
