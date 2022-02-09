package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;
import scratch.UCERF3.inversion.U3InversionTargetMFDs;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

public class SectionMFDPlot extends AbstractPlot {
	
	private RSQSimSubSectionMapper mapper;
	private FaultSystemSolution comparisonSol;
	
	private Map<Integer, IncrementalMagFreqDist> parentSectMFDs;
	private Map<Integer, String> parentSectNames;
	
	private static final double minMag = 5.05;
	private static final int numMag = 41;
	private static final double deltaMag = 0.1;
	
	private EvenlyDiscretizedFunc mfdXVals;
	
	private double minMagEncountered = Double.POSITIVE_INFINITY;

	public SectionMFDPlot(RSQSimSubSectionMapper mapper, FaultSystemSolution comparisonSol) {
		this.mapper = mapper;
		this.comparisonSol = comparisonSol;
		
		mfdXVals = new EvenlyDiscretizedFunc(minMag, numMag, deltaMag);
		
		parentSectMFDs = new HashMap<>();
		parentSectNames = new HashMap<>();
		for (FaultSection sect : mapper.getSubSections()) {
			Integer parentID = sect.getParentSectionId();
			if (!parentSectMFDs.containsKey(parentID)) {
				parentSectMFDs.put(parentID, new IncrementalMagFreqDist(minMag, numMag, deltaMag));
				parentSectNames.put(parentID, sect.getParentSectionName());
			}
		}
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		minMagEncountered = Math.min(mag, minMagEncountered);
		if (mag < minMag-0.5*deltaMag)
			return;
		
		List<List<SubSectionMapping>> bundled = mapper.getFilteredSubSectionMappings((RSQSimEvent)e);
		if (bundled.isEmpty())
			bundled = mapper.getAllSubSectionMappings((RSQSimEvent)e);
		
		HashSet<Integer> parentIDs = new HashSet<>();
		for (List<SubSectionMapping> bundle : bundled)
			for (SubSectionMapping mapping : bundle)
				parentIDs.add(mapping.getSubSect().getParentSectionId());
		
		int mfdIndex = mfdXVals.getClosestXIndex(mag);
		for (Integer parentID : parentIDs) {
			IncrementalMagFreqDist mfd = parentSectMFDs.get(parentID);
			mfd.add(mfdIndex, 1d);
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		int minMagIndex = mfdXVals.getClosestXIndex(minMagEncountered);
		Range xRange = new Range(mfdXVals.getX(minMagIndex)-0.5*deltaMag, mfdXVals.getMaxX()+0.5*deltaMag);
		
		GridSourceProvider gridProv = null;
		FaultGridAssociations polyMgr = null;
		Map<Integer, List<FaultSection>> parentToSubSects = null;
		if (comparisonSol != null && comparisonSol.getGridSourceProvider() != null) {
			polyMgr = FaultPolyMgr.create(comparisonSol.getRupSet().getFaultSectionDataList(), U3InversionTargetMFDs.FAULT_BUFFER);
			gridProv = comparisonSol.getGridSourceProvider();
			parentToSubSects = new HashMap<>();
			for (FaultSection sect : comparisonSol.getRupSet().getFaultSectionDataList()) {
				List<FaultSection> parentSects = parentToSubSects.get(sect.getParentSectionId());
				if (parentSects == null) {
					parentSects = new ArrayList<>();
					parentToSubSects.put(sect.getParentSectionId(), parentSects);
				}
				parentSects.add(sect);
			}
		}
		
		for (Integer parentID : parentSectMFDs.keySet()) {
			IncrementalMagFreqDist incrMFD = parentSectMFDs.get(parentID);
			incrMFD.scale(1d/getCurrentDurationYears()); // annualize
			EvenlyDiscretizedFunc cumMFD = incrMFD.getCumRateDistWithOffset();
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			incrMFD.setName("Incremental");
			funcs.add(incrMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
			cumMFD.setName("Cumulative");
			funcs.add(cumMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, getPrimaryColor()));
			
			double minY = Double.POSITIVE_INFINITY;
			double maxY = 0d;
			for (DiscretizedFunc func : funcs) {
				for (Point2D pt : func) {
					if (pt.getY() > 0) {
						minY = Math.min(minY, pt.getY());
						maxY = Math.max(maxY, pt.getY());
					}
				}
			}
			
			if (comparisonSol != null) {
				IncrementalMagFreqDist compIncrMFD = comparisonSol.calcParticipationMFD_forParentSect(
						parentID, minMag, mfdXVals.getMaxX(), numMag);
				if (gridProv != null) {
					// add gridded
					for (FaultSection sect : parentToSubSects.get(parentID)) {
						Map<Integer, Double> nodeFracts = polyMgr.getNodeFractions(sect.getSectionId());
						for (int nodeIndex : nodeFracts.keySet()) {
							IncrementalMagFreqDist nodeMFD = gridProv.getMFD_SubSeisOnFault(nodeIndex);
							Preconditions.checkState((float)nodeMFD.getDelta() == (float)deltaMag);
							int indexOffset = nodeMFD.getClosestXIndex(minMag);
							Preconditions.checkState((float)nodeMFD.getX(indexOffset) == (float)minMag);
							double fract = nodeFracts.get(nodeIndex);
							for (int i=0; i<numMag; i++) {
								int index = i+indexOffset;
								if (index == nodeMFD.size())
									break;
								compIncrMFD.add(i, nodeMFD.getY(index)*fract);
							}
						}
					}
				}
				EvenlyDiscretizedFunc compMFD = compIncrMFD.getCumRateDistWithOffset();
				
				compMFD.setName("UCERF3 Cumulative");
				funcs.add(compMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, getComparableColor()));
				
				// only increase max with comparison
				maxY = Math.max(maxY, compMFD.getMaxY());
				if (cumMFD.getMaxY() == 0d) {
					// also do min, we don't have anything on this fault
					for (Point2D pt : compMFD)
						if (pt.getY() > 0)
							minY = Math.min(minY, pt.getY());
				}
			}
			
			// buffer them a little
			minY = Math.pow(10, Math.log10(minY)-0.3);
			maxY = Math.pow(10, Math.log10(maxY)+0.3);
			
			Range yRange = calcEncompassingLog10Range(minY, maxY);
			
			String title = getParentSectName(parentID)+" MFDs";
			
			PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", "Rate (1/yr)");
			plot.setLegendVisible(true);
			
			String prefix = getParentSectPrefix(parentID);
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(plot, false, true, xRange, yRange);
			gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
			gp.saveAsPNG(new File(getOutputDir(), prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(getOutputDir(), prefix+".pdf").getAbsolutePath());
			
			CSVFile<String> csv = new CSVFile<>(true);
			if (comparisonSol == null)
				csv.addLine("Magnitude", "Incremental Rate", "Cumulative Rate");
			else
				csv.addLine("Magnitude", "Incremental Rate", "Cumulative Rate", "UCERF3 Cumulative Rate");
			for (int i=minMagIndex; i<numMag; i++) {
				List<String> line = new ArrayList<>();
				line.add((float)mfdXVals.getX(i)+"");
				for (DiscretizedFunc func : funcs)
					line.add((float)func.getY(i)+"");
				csv.addLine(line);
			}
			csv.writeToFile(new File(getOutputDir(), prefix+".csv"));
		}
	}
	
	public String getParentSectName(int parentID) {
		return parentSectNames.get(parentID);
	}
	
	public String getParentSectPrefix(int parentID) {
		String prefix = this.getOutputPrefix();
		if (prefix == null || prefix.isEmpty())
			prefix = "";
		else
			prefix += "_";
		String parentPrefix = getParentSectName(parentID).replaceAll("\\W+", "_");
		while (parentPrefix.contains("__"))
			parentPrefix.replaceAll("__", "_");
		while (parentPrefix.endsWith("_"))
			parentPrefix = parentPrefix.substring(0, parentPrefix.length()-1);
		return prefix+parentPrefix;
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		File gitDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		
		double skipYears = 5000;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_3062.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.JG_tuneBase1m.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance(baseDir);
		
		File catalogOutputDir = new File(gitDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File outputDir = new File(catalogOutputDir, "parent_sect_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		FaultSystemSolution compSol = catalog.getComparisonSolution();
		
		SectionMFDPlot plot = new SectionMFDPlot(catalog.getSubSectMapper(), compSol);
		plot.initialize(catalog.getName(), resourcesDir, "");
		plot.setPlotSize(800, 600);
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();
		
		System.out.println("Building markdown...");
		
		// header
		List<String> lines = new ArrayList<>();
		lines.add("# "+catalog.getName()+" Parent Section MFDs");
		lines.add("");
		lines.add("*Sections participate in a rupture if at least "
				+(float)(catalog.getMinSubSectFractForInclusion()*100d)+" % of one of its subsections (by area) ruptures*");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		Map<String, List<Integer>> namedFaultsMap = new HashMap<>(catalog.getFaultModel().getNamedFaultsMapAlt());
		
		List<String> namedFaults = new ArrayList<>(namedFaultsMap.keySet());
		Collections.sort(namedFaults);
		
		HashSet<Integer> parentsLeft = new HashSet<>(plot.parentSectNames.keySet());
		for (String fault : namedFaults) {
			// sort by name
			List<Integer> parentIDs = namedFaultsMap.get(fault);
			List<String> names = new ArrayList<>();
			for (int parentID : parentIDs) {
				names.add(plot.getParentSectName(parentID));
				parentsLeft.remove(parentID);
			}
			namedFaultsMap.put(fault, ComparablePairing.getSortedData(names, parentIDs));
		}
		
		// now handle other
		List<Integer> otherParentIDs = new ArrayList<>(parentsLeft);
		List<String> otherParentNames = new ArrayList<>();
		for (int parentID : otherParentIDs)
			otherParentNames.add(plot.getParentSectName(parentID));
		String otherName = "Other Faults";
		namedFaults.add(otherName);
		namedFaultsMap.put(otherName, ComparablePairing.getSortedData(otherParentNames, otherParentIDs));
		
		for (String faultName : namedFaults) {
			lines.add("## "+faultName);
			lines.add(topLink); lines.add("");
			
			for (int parentID : namedFaultsMap.get(faultName)) {
				String parentName = plot.getParentSectName(parentID);
				lines.add("### "+parentName);
				lines.add(topLink); lines.add("");
				
				String prefix = plot.getParentSectPrefix(parentID);
				lines.add("![MFD](resources/"+prefix+".png)");
				lines.add("");
				CSVFile<String> csv = CSVFile.readFile(new File(resourcesDir, prefix+".csv"), true);
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				for (String column : csv.getLine(0))
					table.addColumn(column);
				table.finalizeLine();
				for (int row=1; row<csv.getNumRows(); row++) {
					boolean allZeros = true;
					for (int col=1; col<csv.getNumCols(); col++)
						allZeros = allZeros && Float.parseFloat(csv.get(row, col)) == 0f;
					if (allZeros)
						continue;
					table.initNewLine();
					for (String column : csv.getLine(row))
						table.addColumn(column);
					table.finalizeLine();
				}
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);

		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(gitDir);

		System.out.println("DONE");
	}

}
