package scratch.kevin.simulators.multiFault;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.IDPairing;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.utils.FaultSystemIO;

public class RupSetConnectionSearch {
	
	private FaultSystemRupSet rupSet;
	private Map<IDPairing, Double> sectDistanceCache;
	
	private RuptureSurface[] sectSurfs;
	
	public static final double MAX_POSSIBLE_JUMP_DEFAULT = 100d;
	private double maxPossibleJump;
	
	public static final boolean CUMULATIVE_JUMPS_DEFAULT = false;
	// if true, find connections via the smallest cumulative jump distance
	// if false, find connections via the smallest individual jump (possibly across multiple clusters)
	private boolean cumulativeJumps;

	public RupSetConnectionSearch(FaultSystemRupSet rupSet) {
		this(rupSet, new HashMap<>());
	}
	
	public RupSetConnectionSearch(FaultSystemRupSet rupSet, Map<IDPairing, Double> sectDistanceCache) {
		this(rupSet, sectDistanceCache, MAX_POSSIBLE_JUMP_DEFAULT, CUMULATIVE_JUMPS_DEFAULT);
	}
	
	public RupSetConnectionSearch(FaultSystemRupSet rupSet, Map<IDPairing, Double> sectDistanceCache,
			double maxPossibleJump, boolean cumulativeJumps) {
		this.rupSet = rupSet;
		if (sectDistanceCache == null)
			sectDistanceCache = new HashMap<>();
		this.sectDistanceCache = sectDistanceCache;
		this.maxPossibleJump = maxPossibleJump; 
		this.cumulativeJumps = cumulativeJumps;
		
		sectSurfs = new RuptureSurface[rupSet.getNumSections()];
		for (int s=0; s<sectSurfs.length; s++)
			// 1km spacing, don't preserve grid spacing exactly, don't reduce upper depth
			sectSurfs[s] = rupSet.getFaultSectionData(s).getFaultSurface(1d, false, false);
	}
	
	private static IDPairing pair(FaultSection sect1, FaultSection sect2) {
		return pair(sect1.getSectionId(), sect2.getSectionId());
	}
	
	private static IDPairing pair(int id1, int id2) {
		if (id1 < id2)
			return new IDPairing(id1, id2);
		return new IDPairing(id2, id1);
	}
	
	private double calcSectDist(FaultSection sect1, FaultSection sect2) {
		return calcSectDist(sect1.getSectionId(), sect2.getSectionId());
	}
	
	private double calcSectDist(int id1, int id2) {
		if (id1 == id2)
			return 0d;
		IDPairing pair = pair(id1, id2);
		Double minDist = sectDistanceCache.get(pair);
		if (minDist != null)
			return minDist;
		
		RuptureSurface surf1 = sectSurfs[id1];
		RuptureSurface surf2 = sectSurfs[id2];
		
		// if the quick distance is less than this value, calculate a full distance
		double quickDistThreshold = 5d*Math.max(surf1.getAveLength(), surf2.getAveLength());
		minDist = Double.POSITIVE_INFINITY;
		
		for (Location loc : surf1.getPerimeter()) {
			minDist = Math.min(minDist, surf2.getQuickDistance(loc));
			if (minDist < quickDistThreshold)
				break;
		}
		
		if (minDist < quickDistThreshold)
			// do the full calculation
			minDist = surf1.getMinDistance(surf2);
		
		sectDistanceCache.putIfAbsent(pair, minDist);
		
		return minDist;
	}
	
	private List<Cluster> calcClusters(List<FaultSection> sects, final boolean debug) {
		List<Cluster> clusters = new ArrayList<>();
		
		Map<Integer, List<FaultSection>> parentsMap = new HashMap<>();
		
		for (FaultSection sect : sects) {
			Integer parentID = sect.getParentSectionId();
			List<FaultSection> parentSects = parentsMap.get(parentID);
			if (parentSects == null) {
				parentSects = new ArrayList<>();
				parentsMap.put(parentID, parentSects);
			}
			parentSects.add(sect);
		}
		
		for (List<FaultSection> parentSects : parentsMap.values()) {
			// TODO this implementation is specific to U3 subsectioning, might need to be generalized
			// for models that allow subsections down dip
			Collections.sort(parentSects, sectIDcomp);
			
			List<FaultSection> curSects = new ArrayList<>();
			int prevID = -2;
			for (int s=0; s<parentSects.size(); s++) {
				FaultSection sect = parentSects.get(s);
				int id = sect.getSectionId();
				if (!curSects.isEmpty() && id != prevID+1) {
					// new cluster within this section
					clusters.add(new Cluster(curSects));
					curSects = new ArrayList<>();
				}
				curSects.add(sect);
				prevID = id;
			}
			clusters.add(new Cluster(curSects));
		}
		
		// calculate connections
		for (int i=0; i<clusters.size(); i++) {
			Cluster c1 = clusters.get(i);
			if (debug) {
				String parentNane = c1.parentName;
				System.out.println("\tCluster "+i+" with "+c1.sects.size()+" sections on "+parentNane);
			}
			for (int j=i+1; j<clusters.size(); j++) {
				Cluster c2 = clusters.get(j);
				// find the closest connection point between these two clusters
				FaultSection closestFrom = null;
				FaultSection closestTo = null;
				double minDist = maxPossibleJump;
				for (FaultSection s1 : c1.sects) {
					for (FaultSection s2 : c2.sects) {
						double dist = calcSectDist(s1, s2);
						if (dist < minDist) {
							closestFrom = s1;
							closestTo = s2;
							minDist = dist;
						}
					}
				}
				if (closestFrom != null) {
					// we have a possible connection
					c1.addPossibleConnection(closestFrom, c2, minDist);
					c2.addPossibleConnection(closestTo, c1, minDist);
					if (debug)
						System.out.println("\t\tConnection ["+closestFrom.getSectionId()+"=>"
								+closestTo.getSectionId()+"] with R="+(float)minDist+" to "+c2.parentName);
				}
			}
		}
		
		return clusters;
	}
	
	private static final Comparator<FaultSection> sectIDcomp = new Comparator<FaultSection>() {
		@Override
		public int compare(FaultSection o1, FaultSection o2) {
			return Integer.compare(o1.getSectionId(), o2.getSectionId());
		}
	};
	
	private class Cluster {
		// sections that belong to this cluster
		private List<FaultSection> sects;
		// map of <destination cluster, section on this cluster which connects to that one>
		private Map<Cluster, FaultSection> clusterExitPoints;
		private List<Cluster> sortedConnections;
		private List<Double> sortedDists;
		
		private final String parentName;
		
		private Cluster(List<FaultSection> sects) {
			this.sects = sects;
			clusterExitPoints = new HashMap<>();
			parentName = sects.get(0).getParentSectionName();
			int parentID = sects.get(0).getParentSectionId();
			for (int s=1; s<sects.size(); s++) {
				int myParentID = sects.get(s).getParentSectionId();
				Preconditions.checkState(myParentID == parentID,
						"Parent ID mismatch: %s != %s", parentID, myParentID);
			}
			
			sortedConnections = new ArrayList<>();
			sortedDists = new ArrayList<>();
		}
		
		private void addPossibleConnection(FaultSection from, Cluster oCluster, double distance) {
			Preconditions.checkState(!clusterExitPoints.containsKey(oCluster));
			clusterExitPoints.put(oCluster, from);
			int insIndex = Collections.binarySearch(sortedDists, distance);
			if (insIndex < 0) {
				// index = -ins_pt - 1
				// index + 1 = -ins_pt
				// ins_pt = -(index + 1)
				insIndex = -(insIndex + 1);
			}
			sortedConnections.add(insIndex, oCluster);
			sortedDists.add(insIndex, distance);
		}
	}
	
	private class ClusterPath {
		private final Cluster start;
		private final Cluster target;
		
		private final HashSet<Cluster> availableClusters;
		
		private final Cluster[] path;
		private final FaultSection[] exitPoints;
		private final FaultSection[] entryPoints;
		
		private double dist;

		public ClusterPath(Cluster start, Cluster target, HashSet<Cluster> availableClusters) {
			this(start, target, availableClusters,
					new Cluster[] {start}, // path starts with just this cluster
					new FaultSection[0], // exit point not yet known
					new FaultSection[] { null }, // entry point to first cluster is null
					0d); // start at zero dist (no jumps yet)
		}
		
		public ClusterPath(Cluster start, Cluster target, HashSet<Cluster> availableClusters,
				Cluster[] path, FaultSection[] exitPoints, FaultSection[] entryPoints, double dist) {
			this.start = start;
			this.target = target;
			this.availableClusters = new HashSet<>(availableClusters);
			
			this.path = path;
			this.exitPoints = exitPoints;
			this.entryPoints = entryPoints;
			
			this.dist = dist;
		}
		
		public ClusterPath take(Cluster to) {
			Preconditions.checkState(!isComplete());
			
			HashSet<Cluster> newAvailClusters = new HashSet<>(availableClusters);
			Preconditions.checkState(newAvailClusters.remove(to));
			
			Cluster from = path[path.length-1];
			Cluster[] newPath = Arrays.copyOf(path, path.length+1);
			newPath[path.length] = to;
			
			FaultSection exitPoint = from.clusterExitPoints.get(to);
			Preconditions.checkNotNull(exitPoint);
			FaultSection[] newExitPoints = Arrays.copyOf(exitPoints, exitPoints.length+1);
			newExitPoints[exitPoints.length] = exitPoint;
			
			FaultSection entryPoint = to.clusterExitPoints.get(from);
			Preconditions.checkNotNull(entryPoint);
			FaultSection[] newEntryPoints = Arrays.copyOf(entryPoints, entryPoints.length+1);
			newEntryPoints[entryPoints.length] = entryPoint;
			
			double newDist;
			if (cumulativeJumps)
				newDist = dist + calcSectDist(exitPoint, entryPoint);
			else
				newDist = Math.max(dist, calcSectDist(exitPoint, entryPoint));
			
			return new ClusterPath(start, target, newAvailClusters, newPath, newExitPoints, newEntryPoints, newDist);
		}
		
		public boolean isComplete() {
			return path[path.length-1] == target;
		}
		
		public String toString() {
			String str = "";
			double myDist = 0;
			for (int p=0; p<path.length; p++) {
				if (p == 0)
					str += "[";
				else {
					// close out previous
					int exitID = exitPoints[p-1].getSectionId();
					int entryID = entryPoints[p].getSectionId();
					if (cumulativeJumps)
						myDist += calcSectDist(exitPoints[p-1], entryPoints[p]);
					else
						myDist = calcSectDist(exitPoints[p-1], entryPoints[p]);
					str += "; "+exitID+"] => ["+entryID+"; R="+distDF.format(myDist)+"; ";
				}
				str += path[p].parentName;
			}
			str += "]";
			
			if (isComplete())
				str += " COMPLETE R="+distDF.format(dist);
			else
				str += " INCOMPLETE (target: "+target.parentName+") R="+distDF.format(dist);
			return str;
		}
	}
	
	private static final DecimalFormat distDF = new DecimalFormat("0.0");
	
	private class PathResult {
		private ClusterPath shortestPath;
		private int completePathCount = 0;
		
		private void addPath(ClusterPath path, final boolean debug) {
			Preconditions.checkState(path.isComplete());
			if (shortestPath == null || shortestPath.dist > path.dist) {
				if (debug)
					System.out.println("\t\t\t\tNew shortest with R="
							+distDF.format(path.dist)+": "+path);
				shortestPath = path;
			}
			completePathCount++;
		}
	}
	
	private void pathSearch(final ClusterPath basePath, final PathResult result, final boolean debug) {
		Preconditions.checkState(!basePath.isComplete());
		
		Cluster from = basePath.path[basePath.path.length-1];
		
		// search in order of increasing distance
		for (Cluster to : from.sortedConnections) {
			if (!basePath.availableClusters.contains(to))
				continue;
//		for (Cluster to : basePath.availableClusters) {
//			if (!from.clusterExitPoints.containsKey(to))
//				continue;
			ClusterPath path = basePath.take(to);
			if (debug)
				System.out.println("\t\t\t"+path);
//			if (debug) System.out.println("\t\t\t\tTaking "+sect.getSectionId()+" complete="+path.isComplete()
//							+", dist="+(float)path.cumulativeDistance);
			if (path.isComplete()) {
				result.addPath(path, debug);
			} else {
				if (result.shortestPath != null && path.dist >= result.shortestPath.dist)
					// already worse than our current best, stop here
					continue;
				pathSearch(path, result, debug);
			}
		}
	}
	
	public HashSet<IDPairing> calcConnections(int rupIndex) {
		return calcConnections(rupIndex, false);
	}
	
	public HashSet<IDPairing> calcConnections(int rupIndex, final boolean debug) {
		List<FaultSection> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
		
		if (debug) System.out.println("Building clusters for "+rupIndex);
		
		// calculate clusters (between which there may be connections
		List<Cluster> clusters = calcClusters(sects, debug);
		
		if (debug) System.out.println("\tHave "+clusters.size()+" clusters");
		
		HashSet<IDPairing> connections = new HashSet<>();
		
		if (debug) System.out.println("Searching for connections...");
		
		int numCompletePaths = 0;
		
		for (int i=0; i<clusters.size(); i++) {
			Cluster from = clusters.get(i);
			HashSet<Cluster> availableClusters = new HashSet<>(clusters);
			availableClusters.remove(from);
			
			if (debug) System.out.println("\tFrom cluster "+0+", "+from.parentName);
			
			// TODO: add check for no connections possible?
			
			for (int j=i+1; j<clusters.size(); j++) {
				Cluster target = clusters.get(j);
				ClusterPath basePath = new ClusterPath(from, target, availableClusters);
				
				if (debug) System.out.println("\t\tSearching to cluster "+j+", "+target.parentName);
				
				PathResult result = new PathResult();
				pathSearch(basePath, result, debug);
				
				ClusterPath shortest = result.shortestPath;
				if (shortest != null) {
					// we have a valid path
					
					// only keep the first jump (as others are constrained and may not represent shortest path)
					FaultSection sect1 = shortest.exitPoints[0];
					FaultSection sect2 = shortest.entryPoints[1];
					connections.add(pair(sect1, sect2));
					if (debug) System.out.println("\t\t\tShortest path: "+shortest+" (of "+result.completePathCount+")");
				} else if (debug) {
					System.out.println("\t\t\tNo valid path found");
				}
				numCompletePaths += result.completePathCount;
			}
		}
		
		if (debug)
			System.out.println("Found "+connections.size()+" connections (searched "+numCompletePaths+" full paths)");
		
		return connections;
	}
	
	public void plotConnections(File outputDir, String prefix, int rupIndex) throws IOException {
		
	}
	
	public void plotConnections(File outputDir, String prefix, int rupIndex, Set<IDPairing> highlightConn, String highlightName)
			throws IOException {
		HashSet<IDPairing> connections = calcConnections(rupIndex, true);
		List<FaultSection> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
		List<Cluster> clusters = calcClusters(sects, false);
		HashSet<Integer> parentIDs = new HashSet<>();
		for (FaultSection sect : sects)
			parentIDs.add(sect.getParentSectionId());
		
		Color connectedColor = Color.GREEN.darker();
		Color highlightColor = Color.RED.darker();
		Color faultColor = Color.DARK_GRAY;
		Color faultOutlineColor = Color.LIGHT_GRAY;
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Map<Integer, Location> middles = new HashMap<>();
		
		for (int s=0; s<sects.size(); s++) {
			FaultSection sect = sects.get(s);
			RuptureSurface surf = sect.getFaultSurface(1d);
			
			XY_DataSet trace = new DefaultXY_DataSet();
			for (Location loc : surf.getEvenlyDiscritizedUpperEdge())
				trace.set(loc.getLongitude(), loc.getLatitude());
			
			if (sect.getAveDip() != 90d) {
				XY_DataSet outline = new DefaultXY_DataSet();
				LocationList perimeter = surf.getPerimeter();
				for (Location loc : perimeter)
					outline.set(loc.getLongitude(), loc.getLatitude());
				Location first = perimeter.first();
				outline.set(first.getLongitude(), first.getLatitude());
				
				funcs.add(0, outline);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, faultOutlineColor));
			}
			
			middles.put(sect.getSectionId(), GriddedSurfaceUtils.getSurfaceMiddleLoc(surf));
			
			if (s == 0)
				trace.setName("Fault Sections");
			
			funcs.add(trace);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, faultColor));
		}
		
		boolean first = true;
		double maxDist = 0d;
		for (IDPairing connection : connections) {
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			maxDist = Math.max(maxDist, calcSectDist(connection.getID1(), connection.getID2()));
			
			if (first) {
				xy.setName("Connections");
				first = false;
			}
			
			Location loc1 = middles.get(connection.getID1());
			Location loc2 = middles.get(connection.getID2());
			
			xy.set(loc1.getLongitude(), loc1.getLatitude());
			xy.set(loc2.getLongitude(), loc2.getLatitude());
			
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, connectedColor));
		}
		
		if (highlightConn != null) {
			boolean firstHighlight = true;
			for (IDPairing connection : connections) {
				if (!highlightConn.contains(connection))
					continue;
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				if (firstHighlight) {
					xy.setName(highlightName);
					firstHighlight = false;
				}
				
				Location loc1 = middles.get(connection.getID1());
				Location loc2 = middles.get(connection.getID2());
				
				xy.set(loc1.getLongitude(), loc1.getLatitude());
				xy.set(loc2.getLongitude(), loc2.getLatitude());
				
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, highlightColor));
			}
		}
		
		for (XY_DataSet xy : funcs) {
			for (Point2D pt : xy) {
				latTrack.addValue(pt.getY());
				lonTrack.addValue(pt.getX());
			}
		}
		
		XY_DataSet[] outlines = PoliticalBoundariesData.loadCAOutlines();
		PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(PlotLineType.SOLID, (float)1d, Color.GRAY);
		
		for (XY_DataSet outline : outlines) {
			funcs.add(outline);
			chars.add(outlineChar);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Rupture "+rupIndex+" Connections", "Longitude", "Latitude");
		spec.setLegendVisible(true);
		
		Range xRange = new Range(lonTrack.getMin()-0.5, lonTrack.getMax()+0.5);
		Range yRange = new Range(latTrack.getMin()-0.5, latTrack.getMax()+0.5);
		
		double annX = xRange.getLowerBound() + 0.975*xRange.getLength();
		double annYmult = 0.975d;
		double deltaAnnYmult = 0.05;
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 22);
		
		double annY = yRange.getLowerBound() + annYmult*yRange.getLength();
		XYTextAnnotation cAnn = new XYTextAnnotation(
				clusters.size()+" clusters on "+parentIDs.size()+" parent sects",
				annX, annY);
		cAnn.setFont(annFont);
		cAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(cAnn);
		
		annYmult -= deltaAnnYmult;
		annY = yRange.getLowerBound() + annYmult*yRange.getLength();
		XYTextAnnotation jumpAnn = new XYTextAnnotation(
				connections.size()+" connections (max dist: "+distDF.format(maxDist)+")",
				annX, annY);
		jumpAnn.setFont(annFont);
		jumpAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(jumpAnn);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		double len = Math.max(xRange.getLength(), yRange.getLength());
//		double tick = 2d;
		double tick;
		if (len > 6)
			tick = 2d;
		else if (len > 3)
			tick = 1d;
		else if (len > 1)
			tick = 0.5;
		else
			tick = 0.1;
		TickUnits tus = new TickUnits();
		TickUnit tu = new NumberTickUnit(tick);
		tus.add(tu);
		gp.getXAxis().setStandardTickUnits(tus);
		gp.getYAxis().setStandardTickUnits(tus);
		
		
		File file = new File(outputDir, prefix+".png");
		System.out.println("writing "+file.getAbsolutePath());
		double aspectRatio = yRange.getLength() / xRange.getLength();
		gp.getChartPanel().setSize(800, 200 + (int)(600d*aspectRatio));
		gp.saveAsPNG(file.getAbsolutePath());
	}

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File fssFile = new File("/home/kevin/Simulators/catalogs/rundir4983_stitched/fss/"
				+ "rsqsim_sol_m6.5_skip5000_sectArea0.2.zip");
//		int[] plotIndexes = { 1001, 27845, 173243, 193669 };
//		double debugDist = Double.POSITIVE_INFINITY;
		int[] plotIndexes = {  };
		double debugDist = 30d;
		double maxPossibleJumpDist = MAX_POSSIBLE_JUMP_DEFAULT;
		File outputDir = new File("/tmp/rup_conn_rsqsim");
		
//		File fssFile = new File("/home/kevin/workspace/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/"
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
//		int[] plotIndexes = { 25000, 50000, 75000, 100000, 125000, 150000, 175000, 200000, 225000, 238293, 250000 };
//		double debugDist = Double.POSITIVE_INFINITY;
////		int[] plotIndexes = {};
////		double debugDist = 5d;
//		double maxPossibleJumpDist = 15d;
//		File outputDir = new File("/tmp/rup_conn_u3");
		
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(fssFile);
		
		RupSetConnectionSearch search = new RupSetConnectionSearch(rupSet, new HashMap<>(),
				maxPossibleJumpDist, CUMULATIVE_JUMPS_DEFAULT);
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		if (plotIndexes != null && plotIndexes.length > 0) {
			// just plots
			for (int r : plotIndexes)
				search.plotConnections(outputDir, "rup_"+r, r);
		} else {
			HashSet<IDPairing> allConnections = new HashSet<>();
			
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				if (r % 1000 == 0)
					System.out.println("Calculating for rupture "+r+"/"+rupSet.getNumRuptures()
						+" ("+allConnections.size()+" connections found so far)");
				HashSet<IDPairing> rupConnections = search.calcConnections(r);
				boolean debug = false;
				for (IDPairing pair : rupConnections) {
					double dist = search.calcSectDist(pair.getID1(), pair.getID2());
					if (!allConnections.contains(pair) && dist > debugDist) {
						System.out.println("Pairing "+pair+" has dist="+(float)dist);
						debug = true;
					}
				}
				if (debug)
					search.plotConnections(outputDir, "rup_"+r, r);
				allConnections.addAll(rupConnections);
			}
			
			System.out.println("Found "+allConnections.size()+" total connections");
		}
	}

}
