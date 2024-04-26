package scratch.kevin.nshm23.uncertCorrFigures;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;

public class MultiBranchAverageHazardPlotter {

	public static void main(String[] args) throws IOException {
		System.setProperty("java.awt.headless", "true");
		File dir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3");
		LogicTree<?> tree = LogicTree.read(new File(dir, "logic_tree_full_gridded.json"));
		File resultsFile = new File(dir, "results_hazard_full_gridded.zip");
		
		if (args.length == 1)
			tree = tree.sample(Integer.parseInt(args[0]), true);
		
		File plotsDir = new File(dir, "misc_plots");
		Preconditions.checkState(plotsDir.exists() || plotsDir.mkdir());
		File outputDir = new File(plotsDir, "branch_combination_hazard");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		String mapFileName = "map_pga_TWO_IN_50.txt";
		String mapLabel = "PGA, 2% in 50 Year";
		
		GriddedRegion gridReg = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.1, GriddedRegion.ANCHOR_0_0);
		
		List<LogicTreeNode[]> combinations = new ArrayList<>();
		
		combinations.add(new LogicTreeNode[] { SupraSeisBValues.B_0p0 });
		combinations.add(new LogicTreeNode[] { NSHM23_SegmentationModels.NONE });
		combinations.add(new LogicTreeNode[] { SupraSeisBValues.B_0p0, NSHM23_SegmentationModels.NONE });
		combinations.add(new LogicTreeNode[] { SupraSeisBValues.B_1p0 });
		combinations.add(new LogicTreeNode[] { NSHM23_SegmentationModels.CLASSIC });
		combinations.add(new LogicTreeNode[] { SupraSeisBValues.B_1p0, NSHM23_SegmentationModels.CLASSIC });
		
		double[] mean = new double[gridReg.getNodeCount()];
		double totWeight = 0d;
		
		double[] combWeights = new double[combinations.size()];
		double[][] combMeans = new double[combinations.size()][gridReg.getNodeCount()];
		
		ZipFile zip = new ZipFile(resultsFile);
		
		for (int b=0; b<tree.size(); b++) {
			LogicTreeBranch<?> branch = tree.getBranch(b);
			if (b % 1000 == 0 || (tree.size() <= 1000 && b % 100 == 0))
				System.out.println("Loading branch "+b+"/"+tree.size());
			String fName = branch.buildFileName()+"/"+mapFileName;
			
			ZipEntry entry = zip.getEntry(fName);
			Preconditions.checkNotNull(entry, "Entry is null: %s", fName);
			BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
			
			GriddedGeoDataSet xyz = readMapReader(bRead, gridReg);
			
			double weight = tree.getBranchWeight(branch);
			totWeight += weight;
			
			for (int i=0; i<mean.length; i++)
				if (Double.isFinite(xyz.get(i)))
					mean[i] += weight*xyz.get(i);
			
			for (int c=0; c<combinations.size(); c++) {
				boolean match = true;
				for (LogicTreeNode node : combinations.get(c)) {
					if (!branch.hasValue(node)) {
						match = false;
						break;
					}
				}
				if (match) {
					combWeights[c] += weight;
					for (int i=0; i<mean.length; i++)
						if (Double.isFinite(xyz.get(i)))
							combMeans[c][i] += weight*xyz.get(i);
				}
			}
			
			bRead.close();
		}
		
		zip.close();
		
		scale(mean, 1d/totWeight);
		for (int c=0; c<combinations.size(); c++)
			scale(combMeans[c], 1d/combWeights[c]);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		Color outlineColor = new Color(0, 0, 0, 180);
		Color faultColor = new Color(0, 0, 0, 100);
		
		mapMaker.setFaultSections(NSHM23_FaultModels.WUS_FM_v3.getFaultSections());
		
		float outlineWidth = 2f;
		mapMaker.setRegionOutlineChar(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, outlineColor));
		mapMaker.setPoliticalBoundaryChar(new PlotCurveCharacterstics(PlotLineType.SOLID, outlineWidth, outlineColor));
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, faultColor));
		mapMaker.setSectOutlineChar(null);
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		
		for (int c=0; c<combinations.size(); c++) {
			GriddedGeoDataSet pDiff = new GriddedGeoDataSet(gridReg);
			
			for (int i=0; i<pDiff.size(); i++) {
				double ref = mean[i];
				double test = combMeans[c][i];
				if (Double.isFinite(ref) && Double.isFinite(test) && ref > 0d && test > 0d) {
					pDiff.set(i, 100d*(test - ref)/ref);
				} else {
					pDiff.set(i, Double.NaN);
				}
			}
			
			String name = "";
			String prefix = "";
			for (LogicTreeNode node : combinations.get(c)) {
				if (!name.isBlank()) {
					name += ", ";
					prefix += "_";
				}
				name += node.getShortName();
				prefix += node.getFilePrefix();
			}
			System.out.println("Writing for "+name);
			mapMaker.plotXYZData(pDiff, pDiffCPT, name+" / Mean, % Change, "+mapLabel);
			mapMaker.setWriteGeoJSON(false);
			mapMaker.plot(outputDir, prefix, " ");
		}
	}

	private static GriddedGeoDataSet readMapReader(BufferedReader bRead, GriddedRegion gridReg) throws IOException {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		String line = bRead.readLine();
		int index = 0;
		while (line != null) {
			line = line.trim();
			if (!line.startsWith("#")) {
				StringTokenizer tok = new StringTokenizer(line);
				double lon = Double.parseDouble(tok.nextToken());
				double lat = Double.parseDouble(tok.nextToken());
				double val = Double.parseDouble(tok.nextToken());
				Location loc = new Location(lat, lon);
				Preconditions.checkState(LocationUtils.areSimilar(loc, gridReg.getLocation(index)));
				xyz.set(index++, val);
			}
			line = bRead.readLine();
		}
		Preconditions.checkState(index == gridReg.getNodeCount());
		return xyz;
	}
	
	private static void scale(double[] map, double scalar) {
		for (int i=0; i<map.length; i++)
			map[i] *= scalar;
	}
}
