package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SlipRatePlots;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels.MinisectionSlipRecord;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class DefModelsPageGen {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/markdown/nshm23-misc/deformation_models");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# NSHM23 Deformation Models");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		CPT slipCPT = SlipRatePlots.linearSlipCPT(40d);
		CPT logSlipCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 2);
		CPT logRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-5d, 5d);
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v2;
		
		List<? extends FaultSection> fullSects = fm.getFaultSections();
		
		NSHM23_DeformationModels avgDM = NSHM23_DeformationModels.AVERAGE;
		
		NSHM23_DeformationModels[] dms = {
				NSHM23_DeformationModels.GEOLOGIC,
				NSHM23_DeformationModels.EVANS,
				NSHM23_DeformationModels.POLLITZ,
				NSHM23_DeformationModels.SHEN_BIRD,
				NSHM23_DeformationModels.ZENG,
		};
		
		Map<Integer, List<MinisectionSlipRecord>> avgMinis = avgDM.getMinisections(fm);
		List<? extends FaultSection> avgSubSects = avgDM.build(fm);
		List<Map<Integer, List<MinisectionSlipRecord>>> dmMinis = new ArrayList<>();
		List<List<? extends FaultSection>> dmSubSects = new ArrayList<>();
		
		for (NSHM23_DeformationModels dm : dms) {
			dmMinis.add(dm.getMinisections(fm));
			dmSubSects.add(dm.build(fm));
		}
		
		lines.add("## Minisections");
		lines.add(topLink); lines.add("");
		
		lines.add("Slip rates were assigned by the deformation modelers for each fault minisection, defined "
				+ "as straight lenghts of fault between locations on the fault trace of no more than 15 km. "
				+ "This section gives slip rates for each model on those original minisections.");
		lines.add("");
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(NSHM23_RegionLoader.loadFullConterminousWUS());
		mapMaker.setRegionOutlineChar(null);
		mapMaker.setSectOutlineChar(null);
		
		List<FaultSection> minisections = new ArrayList<>();
		List<Double> avgMiniSlips = new ArrayList<>();
		List<List<Double>> dmMiniSlips = new ArrayList<>();
		List<List<FaultSection>> dmMiniSects = new ArrayList<>();
		for (int d=0; d<dms.length; d++) {
			dmMiniSlips.add(new ArrayList<>());
			dmMiniSects.add(new ArrayList<>());
		}
		
		for (FaultSection fullSect : fullSects) {
			List<MinisectionSlipRecord> minis = avgMinis.get(fullSect.getSectionId());
			if (minis == null || minis.isEmpty())
				continue;
			for (int m=0; m<minis.size(); m++) {
				FaultSection sect = miniToSect(minis.get(m), fullSect, minisections.size(), m);
				minisections.add(sect);
				avgMiniSlips.add(minis.get(m).slipRate);
				for (int d=0; d<dms.length; d++) {
					MinisectionSlipRecord dmMini = dmMinis.get(d).get(fullSect.getSectionId()).get(m);
					dmMiniSlips.get(d).add(dmMini.slipRate);
					FaultSection dmSect = sect.clone();
					dmSect.setAveSlipRate(dmMini.slipRate);
					dmSect.setSlipRateStdDev(dmMini.slipRateStdDev);
					dmSect.setAveRake(dmMini.rake);
					dmMiniSects.get(d).add(dmSect);
				}
			}
		}
		
		File csvFile = new File(resourcesDir, "minisect_slips.csv");
//		CSVFile<String> csv = buildSlipsCSV(minisections, avgMiniSlips, dmMiniSlips, dms, true);
//		csv.writeToFile(csvFile);
		NSHM23_DeformationModels.writeMinisectionCSV(csvFile);
//		lines.add("__Download CSV File: ["+csvFile.getName()+"]("+resourcesDir.getName()+"/"+csvFile.getName()+")__");
//		lines.add("");
		
		TableBuilder downloadsTable = buildDownloadsTable(minisections, dmMiniSects, dms, csvFile, resourcesDir, "minis");
		lines.add("__Downloads__");
		lines.add("");
		lines.addAll(downloadsTable.build());
		lines.add("");
		lines.add("__Plots__");
		lines.add("");
		
		TableBuilder table = buildSlipsTable(minisections, avgMiniSlips, dmMiniSlips, dms, mapMaker,
				slipCPT, logSlipCPT, diffCPT, logRatioCPT, resourcesDir, "minis");
		
		lines.addAll(table.build());
		
		lines.add("## Subsection Mapped");
		lines.add(topLink); lines.add("");
		
		lines.add("For use in the inversion, faults were divided into equal length subsections, each with "
				+ "length approximately equal to half of the down-dip width. This section gives slip rates for "
				+ "each model mapped onto those subsections. Values are present both before and after slip "
				+ "rate reductions due to creep (which occurs when the creep rate exceeds "
				+ new DecimalFormat("0.#%").format(NSHM23_DeformationModels.ASEIS_CEILING)+" of the slip rate).");
		lines.add("");
		
		for (boolean creepReduced : new boolean[] {false, true}) {
			String prefix;
			if (creepReduced) {
				lines.add("### Subsection Mapped, Creep Reduced");
				prefix = "subsects_creep_reduced";
			} else {
				lines.add("### Subsection Mapped, Original");
				prefix = "subsects_original";
			}
			lines.add(topLink); lines.add("");
			
			List<Double> avgSlips = new ArrayList<>();
			List<List<Double>> dmSlips = new ArrayList<>();
			for (int d=0; d<dms.length; d++)
				dmSlips.add(new ArrayList<>());
			for (int s=0; s<avgSubSects.size(); s++) {
				FaultSection avgSect = avgSubSects.get(s);
				avgSlips.add(creepReduced ? avgSect.getReducedAveSlipRate() : avgSect.getOrigAveSlipRate());
				for (int d=0; d<dms.length; d++) {
					FaultSection dmSect = dmSubSects.get(d).get(s);
					dmSlips.get(d).add(creepReduced ? dmSect.getReducedAveSlipRate() : dmSect.getOrigAveSlipRate());
				}
			}
			
			CSVFile<String> csv = buildSlipsCSV(avgSubSects, avgSlips, dmSlips, dms, false);
			
			csvFile = new File(resourcesDir, prefix+".csv");
			csv.writeToFile(csvFile);
//			lines.add("__Download CSV File: ["+csvFile.getName()+"]("+resourcesDir.getName()+"/"+csvFile.getName()+")__");
//			lines.add("");
			
			downloadsTable = buildDownloadsTable(creepReduced ? null : stripGenerics(avgSubSects),
					creepReduced ? null : stripGenerics2(dmSubSects), dms, csvFile, resourcesDir, prefix);
			lines.add("__Downloads__");
			lines.add("");
			lines.addAll(downloadsTable.build());
			lines.add("");
			lines.add("__Plots__");
			lines.add("");
			
			table = buildSlipsTable(avgSubSects, avgSlips, dmSlips, dms, mapMaker,
					slipCPT, logSlipCPT, diffCPT, logRatioCPT, resourcesDir, prefix);
			
			lines.addAll(table.build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static FaultSection miniToSect(MinisectionSlipRecord mini, FaultSection refSect, int index, int minisectionIndex) {
		FaultSectionPrefData sect = new FaultSectionPrefData();
		sect.setSectionId(index);
		sect.setSectionName(refSect.getSectionName()+" Minisection "+minisectionIndex);
		sect.setParentSectionId(refSect.getSectionId());
		sect.setParentSectionName(refSect.getSectionName());
		
		FaultTrace trace = new FaultTrace(sect.getSectionName());
		trace.add(mini.startLoc);
		trace.add(mini.endLoc);
		sect.setFaultTrace(trace);
		
		sect.setAveUpperDepth(refSect.getOrigAveUpperDepth());
		sect.setAveLowerDepth(refSect.getAveLowerDepth());
		sect.setAveDip(refSect.getAveDip());
		sect.setDipDirection(refSect.getDipDirection());
		sect.setAveRake(refSect.getAveRake());
		sect.setAveSlipRate(mini.slipRate);
		sect.setSlipRateStdDev(mini.slipRateStdDev);
		
		return sect;
	}
	
	private static TableBuilder buildSlipsTable(List<? extends FaultSection> sects,
			List<Double> avgSlipRates, List<List<Double>> dmSlipRates, NSHM23_DeformationModels[] dms,
			GeographicMapMaker mapMaker, CPT slipCPT, CPT logSlipCPT, CPT diffCPT, CPT logRatioCPT,
			File resourcesDir, String prefix) throws IOException {
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		mapMaker.clearFaultSections();
		mapMaker.setFaultSections(sects);
		mapMaker.setWriteGeoJSON(true);
		
		for (int d=-1; d<dms.length; d++) {
			List<Double> rates = d < 0 ? avgSlipRates : dmSlipRates.get(d);
			String name = d < 0 ? "Average" : dms[d].getShortName();
			String dmPrefix = prefix+"_"+(d < 0 ? "avg" : dms[d].getFilePrefix());
			
			List<Double> logRates = new ArrayList<>(rates.size());
			for (double rate : rates)
				logRates.add(Math.log10(rate));
			
			mapMaker.plotSectScalars(rates, slipCPT, name+" Slip Rates (mm/yr)");
			mapMaker.plot(resourcesDir, dmPrefix, " ");
			mapMaker.plotSectScalars(logRates, logSlipCPT, "Log10("+name+" Slip Rates) (mm/yr)");
			mapMaker.plot(resourcesDir, dmPrefix+"_log", " ");
			
			table.addLine(MarkdownUtils.boldCentered(name), "");
			table.initNewLine();
			table.addColumn("!["+name+"]("+resourcesDir.getName()+"/"+dmPrefix+".png)");
			table.addColumn("!["+name+"]("+resourcesDir.getName()+"/"+dmPrefix+"_log.png)");
			table.finalizeLine().initNewLine();
			String relPath = resourcesDir.getName();
			table.addColumn(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON", relPath+"/"+dmPrefix+".geojson"));
//					+" "+"[Download GeoJSON]("+relPath+"/"+dmPrefix+".geojson)");
			table.addColumn(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON", relPath+"/"+dmPrefix+"_log.geojson"));
//					+" "+"[Download GeoJSON]("+relPath+"/"+dmPrefix+"_log.geojson)");
			table.finalizeLine();
			if (d >= 0) {
				table.initNewLine();
				List<Double> diffs = new ArrayList<>();
				List<Double> ratios = new ArrayList<>();
				for (int i=0; i<sects.size(); i++) {
					double rate = rates.get(i);
					double avgRate = avgSlipRates.get(i);
					diffs.add(rate - avgRate);
					ratios.add(Math.log10(rate/avgRate));
				}
				mapMaker.plotSectScalars(diffs, diffCPT, name+" - Average (mm/yr)");
				mapMaker.plot(resourcesDir, dmPrefix+"_diff", " ");
				mapMaker.plotSectScalars(ratios, logRatioCPT, "Log10("+name+"/Average)");
				mapMaker.plot(resourcesDir, dmPrefix+"_ratio", " ");
				
				table.addColumn("!["+name+"]("+resourcesDir.getName()+"/"+dmPrefix+"_diff.png)");
				table.addColumn("!["+name+"]("+resourcesDir.getName()+"/"+dmPrefix+"_ratio.png)");
				table.finalizeLine().initNewLine();
				table.addColumn(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON", relPath+"/"+dmPrefix+"_diff.geojson"));
//						+" "+"[Download GeoJSON]("+relPath+"/"+dmPrefix+"_diff.geojson)");
				table.addColumn(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON", relPath+"/"+dmPrefix+"_ratio.geojson"));
//						+" "+"[Download GeoJSON]("+relPath+"/"+dmPrefix+"_ratio.geojson)");
				table.finalizeLine();
			}
		}
		
		return table;
	}
	
	private static CSVFile<String> buildSlipsCSV(List<? extends FaultSection> sects,
			List<Double> avgSlipRates, List<List<Double>> dmSlipRates,
			NSHM23_DeformationModels[] dms, boolean minis) {
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		
		header.add("Section ID");
		header.add("Section Name");
		if (minis)
			header.add("Minisection Number");
		else
			header.add("Subsection Number");
		header.add("Average Slip Rate (mm/yr)");
		for (NSHM23_DeformationModels dm : dms)
			header.add(dm.getShortName()+" Slip Rate (mm/yr)");
		csv.addLine(header);
		
		int curParentID = -1;
		int parentIDStartIndex = -1;
		for (int i=0; i<sects.size(); i++) {
			FaultSection sect = sects.get(i);
			int parentID = sect.getParentSectionId();
			Preconditions.checkState(parentID >= 0);
			Preconditions.checkState(sect.getSectionId() == i);
			if (parentID != curParentID) {
				curParentID = parentID;
				parentIDStartIndex = i;
			}
			int subID = i - parentIDStartIndex;
			List<String> line = new ArrayList<>(header.size());
			line.add(sect.getParentSectionId()+"");
			line.add(sect.getParentSectionName());
			line.add(subID+"");
			line.add(avgSlipRates.get(i).floatValue()+"");
			for (List<Double> rates : dmSlipRates)
				line.add(rates.get(i).floatValue()+"");
			csv.addLine(line);
		}
		
		return csv;
	}
	
	private static List<FaultSection> stripGenerics(List<? extends FaultSection> sects) {
		List<FaultSection> ret = new ArrayList<>();
		for (FaultSection sect : sects)
			ret.add(sect);
		return ret;
	}
	
	private static List<List<FaultSection>> stripGenerics2(List<List<? extends FaultSection>> sects2) {
		List<List<FaultSection>> ret = new ArrayList<>();
		for (List<? extends FaultSection> sects : sects2)
			ret.add(stripGenerics(sects));
		return ret;
	}
	
	private static TableBuilder buildDownloadsTable(List<FaultSection> avgSects,
			List<List<FaultSection>> dmSects, NSHM23_DeformationModels[] dms,
			File csvFile, File resourcesDir, String prefix) throws IOException {
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("Description", "Link");
		table.addLine("Slip Rates CSV File", "["+csvFile.getName()+"]("+resourcesDir.getName()+"/"+csvFile.getName()+")");
		
		if (avgSects != null) {
			File avgSectsFile = new File(resourcesDir, prefix+"_average_slips.geojson");
			GeoJSONFaultReader.writeFaultSections(avgSectsFile, avgSects);
			table.addLine("Average GeoJSON File", "["+avgSectsFile.getName()+"]("+resourcesDir.getName()+"/"+avgSectsFile.getName()+")");
			for (int d=0; d<dms.length; d++) {
				File dmFile = new File(resourcesDir, prefix+"_"+dms[d].getFilePrefix()+"_slips.geojson");
				GeoJSONFaultReader.writeFaultSections(dmFile, dmSects.get(d));
				table.addLine(dms[d].getShortName()+" GeoJSON File", "["+dmFile.getName()+"]("+resourcesDir.getName()+"/"+dmFile.getName()+")");
			}
		}
		
		return table;
	}

}
