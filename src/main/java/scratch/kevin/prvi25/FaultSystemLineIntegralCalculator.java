package scratch.kevin.prvi25;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.title.Title;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainBoundedDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionConnectionList;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;
import com.itextpdf.text.pdf.parser.Line;

import mpi.MaxFloat;

public class FaultSystemLineIntegralCalculator {
	
	private static boolean D = false;

	private List<? extends FaultSection> sects;
	private List<Vector3D> slipVectors;
	private List<FaultTrace> discrUpperEdges;
	private List<Double> slipRates;
	
	public enum VectorComponent {
		FULL_HORIZONTAL("|Full Horizontal|", PlotLineType.SOLID) {
			@Override
			public Vector3D get(Vector3D fullVect, double lineAzimuth) {
				return fullVect;
			}

			@Override
			public double getMagnitude(Vector3D fullVect, double lineAzimuth) {
				return vectorHorzLen(fullVect);
			}
		},
		PARALLEL("Parallel", PlotLineType.DASHED) {
			@Override
			public Vector3D get(Vector3D fullVect, double lineAzimuth) {
				Vector3D lineDirection = directionVector(lineAzimuth);
				Preconditions.checkState((float)vectorHorzLen(lineDirection) == 1f);
				// Parallel component (projection of vector onto the line direction)
				return lineDirection.scalarMultiply(fullVect.dotProduct(lineDirection));
			}

			@Override
			public double getMagnitude(Vector3D fullVect, double lineAzimuth) {
				Vector3D lineDirection = directionVector(lineAzimuth+90);
				Vector3D parallel = get(fullVect, lineAzimuth);
				double cross = lineDirection.getX() * fullVect.getY() - lineDirection.getY() * fullVect.getX();
				double length = vectorHorzLen(parallel);
				if (cross <= 0)
					length = -length;
				return length;
			}
		},
		PERPENDICULAR("Perpendicular", PlotLineType.DOTTED) {
			@Override
			public Vector3D get(Vector3D fullVect, double lineAzimuth) {
				// Perpendicular component (subtract parallel component from the original vector)
				return fullVect.subtract(PARALLEL.get(fullVect, lineAzimuth));
			}

			@Override
			public double getMagnitude(Vector3D fullVect, double lineAzimuth) {
				Vector3D lineDirection = directionVector(lineAzimuth);
				Vector3D perpendicular = get(fullVect, lineAzimuth);
				double cross = lineDirection.getX() * fullVect.getY() - lineDirection.getY() * fullVect.getX();
				double length = vectorHorzLen(perpendicular);
				if (cross >= 0)
					length = -length;
				return length;
			}
		};
		
		public final String label;
		private PlotLineType lineType;

		private VectorComponent(String label, PlotLineType lineType) {
			this.label = label;
			this.lineType = lineType;
		}
		
		public abstract Vector3D get(Vector3D fullVect, double lineAzimuth);
		
		public abstract double getMagnitude(Vector3D fullVect, double lineAzimuth);
	}

	public FaultSystemLineIntegralCalculator(FaultSystemRupSet rupSet) {
		List<Double> slipRates = new ArrayList<>(rupSet.getNumSections());
		for (int s=0; s<rupSet.getNumSections(); s++)
			slipRates.add(rupSet.getSlipRateForSection(s)*1e3); // m/yr -> mm/yr
		init(rupSet.getFaultSectionDataList(), slipRates);
	}
	
	public FaultSystemLineIntegralCalculator(List<? extends FaultSection> sects, boolean reducedSlipRates) {
		List<Double> slips = new ArrayList<>(sects.size());
		for (FaultSection sect : sects)
			slips.add(reducedSlipRates ? sect.getReducedAveSlipRate() : sect.getOrigAveSlipRate());
		init(sects, slips);
	}
	
	public FaultSystemLineIntegralCalculator(List<? extends FaultSection> sects, List<Double> slipRates) {
		init(sects, slipRates);
	}

	private void init(List<? extends FaultSection> sects, List<Double> slipRates) {
		this.sects = sects;
		this.slipRates = slipRates;
		
		discrUpperEdges = new ArrayList<>(sects.size());
		slipVectors = new ArrayList<>(sects.size());
		for (FaultSection sect : sects) {
			RuptureSurface surf = sect.getFaultSurface(1d, false, true);
			discrUpperEdges.add(surf.getEvenlyDiscritizedUpperEdge());
			slipVectors.add(calcHangingWallSlipVector(sect));
		}
	}
	
	private double[] findSectIntersections(Location p0, Location p1) {
		double[] intersections = new double[sects.size()];
		
		double lineLen = LocationUtils.linearDistanceFast(p0, p1);
		
		for (int s=0; s<intersections.length; s++) {
			FaultTrace discrTrace = discrUpperEdges.get(s);
			
			double dist0 = LocationUtils.distanceToLineFast(p0, p1, discrTrace.first());
			double dist1 = LocationUtils.distanceToLineFast(p0, p1, discrTrace.last());
			
			if ((dist0 > 0) != (dist1 > 0)) {
				// the trace crosses this line
				// now find out where
				
				double[] dists = new double[discrTrace.size()];
				double minDist = Double.POSITIVE_INFINITY;
				int minIndex = -1;
				for (int i=0; i<dists.length; i++) {
					dists[i] = LocationUtils.distanceToLineFast(p0, p1, discrTrace.get(i));
					if (Math.abs(dists[i]) < minDist) {
						minDist = Math.abs(dists[i]);
						minIndex = i;
					}
				}
				
				
				
				// we now know the closest location, and it's distance to the line
				Location closestPt = discrTrace.get(minIndex);
				// move it toward the line
				// true if the closest point is to the right of the line
				double lineAzRad = LocationUtils.azimuthRad(p0, p1);
				Location intersectionPt;
				if (dists[minIndex] >= 0) {
					// positive distance means we're on the right side of the line
					// thus we need to move to the left to intersect the line
					double leftAz = lineAzRad - Math.PI*0.5;
					intersectionPt = LocationUtils.location(closestPt, leftAz, minDist); // minDist is absolute already
				} else {
					// negative distance means we're on the left side of the line
					// thus we need to move to the right to intersect the line
					double rightAz = lineAzRad + Math.PI*0.5;
					intersectionPt = LocationUtils.location(closestPt, rightAz, minDist); // minDist is absolute already
				}
				
				double intersectDistPt0 = LocationUtils.horzDistanceFast(p0, intersectionPt);
				double intersectDistPt1 = LocationUtils.horzDistanceFast(p1, intersectionPt);
				if (intersectDistPt0 < intersectDistPt1 && intersectDistPt1 > lineLen)
					// we're before the first point
					intersections[s] = -intersectDistPt0;
				else
					intersections[s] = intersectDistPt0;
			} else {
				intersections[s] = Double.NaN;
			}
		}
		return intersections;
	}
	
	public LineIntegralResult calcLineIntegral(Location startLoc, Location endLoc) {
		if (D) System.out.println("Calculating line intergral from "+startLoc+" to "+endLoc);
		if (D) System.out.println("\tfinding intersecting sections");
		double[] intersections = findSectIntersections(startLoc, endLoc);
		List<Integer> intersectingIndexes = new ArrayList<>();
		List<Double> intersectingDistances = new ArrayList<>();
		for (int s=0; s<intersections.length; s++) {
			if (Double.isFinite(intersections[s])) {
				intersectingIndexes.add(s);
				intersectingDistances.add(intersections[s]);
			}
		}

		double lineLen = LocationUtils.linearDistanceFast(startLoc, endLoc);
		double lineAzRad = LocationUtils.azimuthRad(startLoc, endLoc);
		double lineAz = Math.toDegrees(lineAzRad);
		
		if (intersectingIndexes.isEmpty()) {
			if (D) System.out.println("No intersections found");
			return new LineIntegralResult(startLoc, endLoc, lineLen, lineAz, new Location[0], new double[0], new Vector3D[0]);
		}
		if (D) System.out.println("Found "+intersectingIndexes.size()+" intersections");
		List<Integer> sortedIndexes = ComparablePairing.getSortedData(intersectingDistances, intersectingIndexes);
		
		List<Double> intersectionsWithinLine = new ArrayList<>(sortedIndexes.size());
		for (int index : sortedIndexes) {
			FaultSection sect = sects.get(index);
			double dist = intersections[index];
			if (D) System.out.println("Section "+sect.getName()+" intersects at "+(float)dist+" along transect");
			if (dist > 0d && dist < lineLen) {
				intersectionsWithinLine.add(dist);
			}
		}
		
		List<Double> distsAlongToCalc = new ArrayList<>(2+sortedIndexes.size());
		distsAlongToCalc.add(0d);
		for (int i=1; i<intersectionsWithinLine.size(); i++) {
			double int1 = intersectionsWithinLine.get(i-1);
			double int2 = intersectionsWithinLine.get(i);
			distsAlongToCalc.add(0.5*(int1 + int2));
		}
		distsAlongToCalc.add(lineLen);
		
		Location[] calcLocs = new Location[distsAlongToCalc.size()];
		double[] vectorDistsAlong = new double[calcLocs.length];
		Vector3D[] vectors = new Vector3D[calcLocs.length];
		
		for (int d=0; d<calcLocs.length; d++) {
			double distAlong = distsAlongToCalc.get(d);
			vectorDistsAlong[d] = distAlong;
			Location loc = LocationUtils.location(startLoc, lineAzRad, distAlong);
			if (D) System.out.println("Calculating vector at d="+(float)distAlong+": "+loc);
			Vector3D vectSum = null;
			for (int index : sortedIndexes) {
				double faultDistAlong = intersections[index];
				
				double sectStrike = discrUpperEdges.get(index).getAveStrike();
				Location intersectionLoc = LocationUtils.location(startLoc, lineAzRad, faultDistAlong);
				Location alongStrikeLoc = LocationUtils.location(intersectionLoc, Math.toRadians(sectStrike), 10d);
				boolean hangingWall = LocationUtils.distanceToLineFast(intersectionLoc, alongStrikeLoc, loc) >= 0;
				Vector3D vect = slipVectors.get(index);
				if (!hangingWall)
					// reverse it
					vect = vect.negate();
				// scale by slip rate
				vect = vect.scalarMultiply(slipRates.get(index));
				if (D) System.out.println("\t"+sects.get(index).getName()+": slip vector on my side is "+vect
						+";\tlen="+(float)vectorLen(vect)+";\thorzLen="+(float)vectorHorzLen(vect)+"; az="+(float)vectorAzimuth(vect, loc));
				if (vectSum == null)
					vectSum = vect;
				else
					vectSum = vectSum.add(vect);
			}
			
			if (D) System.out.println("\tFinal slip vector: "+vectSum
					+";\tlen="+(float)vectorLen(vectSum)+";\thorzLen="+(float)vectorHorzLen(vectSum)+"; az="+(float)vectorAzimuth(vectSum, loc));
			
			calcLocs[d] = loc;
			vectors[d] = vectSum;
		}
		
		return new LineIntegralResult(startLoc, endLoc, lineLen, lineAz, calcLocs, vectorDistsAlong, vectors);
	}
	
	public static class LineIntegralResult {
		public final Location startLoc;
		public final Location endLoc;
		public final double length;
		public final double azimuth;
		
		public final Location[] vectorLocs;
		public final double[] vectorDistsAlong;
		public final Vector3D[] vectors;
		
		private LineIntegralResult(Location startLoc, Location endLoc, double length, double azimuth,
				Location[] vectorLocs, double[] vectorDistsAlong, Vector3D[] vectors) {
			super();
			this.startLoc = startLoc;
			this.endLoc = endLoc;
			this.length = length;
			this.azimuth = azimuth;
			this.vectorLocs = vectorLocs;
			this.vectorDistsAlong = vectorDistsAlong;
			this.vectors = vectors;
		}
	}
	
	public static Vector3D calcHangingWallSlipVector(FaultSection sect) {
		double strike = sect.getFaultTrace().getAveStrike();
		FaultTrace lowerTrace = sect.getLowerFaultTrace();
		if (lowerTrace != null) {
			double lowerStrike = lowerTrace.getAveStrike();
			strike = FaultUtils.getAngleAverage(List.of(strike, lowerStrike));
		}
		double dip = sect.getAveDip();
		double rake = sect.getAveRake();
		
		double[] vector = FaultUtils.getSlipVector(new double[] {strike, dip, rake});
		return new Vector3D(vector);
	}
	
	public static double vectorLen(Vector3D vect) {
		return vect.getNorm();
	}
	
	public static double vectorHorzLen(Vector3D vect) {
		return Math.sqrt(vect.getX()*vect.getX() + vect.getY()*vect.getY());
	}
	
	public static double vectorAzimuth(Vector3D vect, Location loc) {
		Location loc2 = LocationUtils.location(loc, Math.PI*0.5, vect.getX());
		loc2 = LocationUtils.location(loc2, 0d, vect.getY());
		return LocationUtils.azimuth(loc, loc2);
	}
	
	public void plotMap(File outputDir, String prefix, String title, boolean labelMagnitudes, boolean plotVectors, double maxSlip,
			List<LineIntegralResult> integrals) throws IOException {
		GeographicMapMaker mapMaker = buildMapPlot(labelMagnitudes, plotVectors, maxSlip, integrals);
		
		mapMaker.plot(outputDir, prefix, title);
	}
	
	public GeographicMapMaker buildMapPlot(boolean labelMagnitudes, boolean plotVectors, double maxSlip,
			List<LineIntegralResult> integrals) throws IOException {
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (FaultSection sect : sects) {
			for (Location loc : sect.getFaultTrace()) {
				latTrack.addValue(loc.lat);
				lonTrack.addValue(loc.lon);
			}
		}
		
		for (LineIntegralResult integral : integrals) {
			latTrack.addValue(integral.startLoc.lat);
			lonTrack.addValue(integral.startLoc.lon);
			latTrack.addValue(integral.endLoc.lat);
			lonTrack.addValue(integral.endLoc.lon);
		}
		
		double annBufferKM = 30d;
		if (labelMagnitudes) {
			double initialLengthScale = Math.max(
					LocationUtils.horzDistanceFast(new Location(latTrack.getCenter(), lonTrack.getMin()),
							new Location(latTrack.getCenter(), lonTrack.getMax())),
					LocationUtils.horzDistanceFast(new Location(latTrack.getMin(), lonTrack.getCenter()),
							new Location(latTrack.getMax(), lonTrack.getCenter())));
			annBufferKM = Math.max(annBufferKM, 0.15*initialLengthScale);
			for (LineIntegralResult integral : integrals) {
				Location annLoc = LocationUtils.location(integral.endLoc, Math.toRadians(integral.azimuth), annBufferKM);
				latTrack.addValue(annLoc.lat);
				lonTrack.addValue(annLoc.lon);
			}
			if (D) System.out.println("Initial length scale: "+(float)initialLengthScale+", annBufferKM: "+(float)annBufferKM);
		}

		double degreeBuffer = 0.5;
		Font numberAnnFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
		Font textAnnFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
		
		Location lowerLeft = new Location(latTrack.getMin()-degreeBuffer, lonTrack.getMin()-degreeBuffer);
		Location upperRight = new Location(latTrack.getMax()+degreeBuffer, lonTrack.getMax()+degreeBuffer);
		
		double middleLat = 0.5*(lowerLeft.lat + upperRight.lat);
		double middleLon = 0.5*(lowerLeft.lon + upperRight.lon);
		
		Location left = new Location(middleLat, lowerLeft.getLongitude());
		Location right = new Location(middleLat, upperRight.getLongitude());
		Location top = new Location(upperRight.getLatitude(), middleLon);
		Location bottom = new Location(lowerLeft.getLatitude(), middleLon);
		double height = LocationUtils.horzDistance(top, bottom);
		double width = LocationUtils.horzDistance(left, right);
		Location centerLoc = new Location(middleLat, middleLon);

		double[] annBufferDists = {
				annBufferKM*1,
				annBufferKM*0.8,
				annBufferKM*0.6,
				annBufferKM*0.4
		};
		
		double maxVectDist = 0.25 * Math.min(height, width);
		for (int i=0; i<integrals.size(); i++) {
			LineIntegralResult integral1 = integrals.get(i);
			for (int j=i+1; j<integrals.size(); j++) {
				LineIntegralResult integral2 = integrals.get(j);
				double interDist = LocationUtils.horzDistanceFast(integral1.startLoc, integral2.startLoc);
				interDist = Math.min(interDist, LocationUtils.horzDistanceFast(integral1.startLoc, integral2.endLoc));
				interDist = Math.min(interDist, LocationUtils.horzDistanceFast(integral1.endLoc, integral2.startLoc));
				interDist = Math.min(interDist, LocationUtils.horzDistanceFast(integral1.endLoc, integral2.endLoc));
				maxVectDist = Math.min(maxVectDist, 0.8*interDist);
			}
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(new Region(lowerLeft, upperRight));
		mapMaker.setFaultSections(sects);
		mapMaker.setWriteGeoJSON(false);
		
//		double maxSlip = slipRates.stream().mapToDouble(S->S.doubleValue()).max().getAsDouble();
//		for (LineIntegralResult result : integrals)
//			for (Vector3D vect : result.vectors)
//				maxSlip = Math.max(maxSlip, vectorHorzLen(vect));
		
//		CPT slipCPT = SlipRatePlots.linearSlipCPT(maxSlip);
		CPT slipCPT;
		if (Double.isFinite(maxSlip) && maxSlip > 0d) {
			slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, maxSlip);
			mapMaker.plotSectScalars(slipRates, slipCPT, plotVectors ? "Slip or Vector Rate (mm/yr)" : "Slip Rate (mm/yr)");
			mapMaker.setScalarThickness(2f);
		} else {
			Preconditions.checkState(!plotVectors, "Can't plot vectors without giving CPT scale");
			slipCPT = null;
		}
		
		List<LocationList> lines = new ArrayList<>();
		List<PlotCurveCharacterstics> lineChars = new ArrayList<>();
		
		List<LocationList> arrows = new ArrayList<>();
		List<Color> arrowColors = new ArrayList<>();
		
		DecimalFormat slipDF;
		if (maxSlip > 5d)
			slipDF = new DecimalFormat("0.0");
		else
			slipDF = new DecimalFormat("0.00");
		DecimalFormat angleDF = new DecimalFormat("0");
		
		LocationList vectorLocs = new LocationList();
		for (LineIntegralResult result : integrals) {
			LocationList line = LocationList.of(result.startLoc, result.endLoc);
			
			if (plotVectors) {
				lines.add(line);
				lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
			} else {
				// plot the lines as direction arrows
				arrows.add(line);
				arrowColors.add(Color.BLACK);
				
				vectorLocs.add(result.startLoc);
				
				if (result == integrals.get(0)) {
					double integralMiddleLat = 0.5*(result.startLoc.lat + result.endLoc.lat);
					double integralMiddleLon = 0.5*(result.startLoc.lon + result.endLoc.lon);
					double furthestLeftLon = Math.max(lowerLeft.lon+0.25, lowerLeft.lon + 0.05*(upperRight.lon - lowerLeft.lon));
					System.out.println("Furthest left="+furthestLeftLon+", integral left="+integralMiddleLon);
					XYTextAnnotation ann = new XYTextAnnotation("Summation Direction", integralMiddleLon, integralMiddleLat);
					ann.setFont(textAnnFont);
					ann.setTextAnchor(TextAnchor.BOTTOM_CENTER);
					ann.setRotationAnchor(TextAnchor.BOTTOM_CENTER);
					double vectorAngle = result.azimuth;
					while (vectorAngle > 180d)
						vectorAngle -= 360d;
					if (vectorAngle > 45d) {
						// pointing to the right enough that we want it on the other side
						vectorAngle -= 180;
						ann.setTextAnchor(TextAnchor.BOTTOM_CENTER);
						ann.setRotationAnchor(TextAnchor.BOTTOM_CENTER);
					}
//					} else if (integralMiddleLon < furthestLeftLon) {
////					} else {
//						// too far to the left, put it on the other side
////						vectorAngle -= 180;
//						ann.setTextAnchor(TextAnchor.BOTTOM_CENTER);
//						ann.setRotationAnchor(TextAnchor.BOTTOM_CENTER);
//					}
					double angleRad = Math.toRadians(vectorAngle + 90);
					ann.setRotationAngle(mapMaker.getRotationAngleCorrectedForAspectRatio(angleRad));
					mapMaker.addAnnotation(ann);
				}
			}
			
			double lineAz = result.azimuth;
			
			for (int i=0; i<result.vectorLocs.length; i++) {
				Location loc = result.vectorLocs[i];
				Vector3D vector = result.vectors[i];
				
				double rawLen = vectorHorzLen(vector);
				double scaledLen = (rawLen/maxSlip) * maxVectDist;
				double az = vectorAzimuth(vector, loc);
				
				if (plotVectors)
					vectorLocs.add(loc);
				
				if (plotVectors) {
					Location loc2 = locConstPlotDist(loc, centerLoc, az, scaledLen);
					LocationList vectLine = LocationList.of(loc, loc2);
					arrows.add(vectLine);
					arrowColors.add(slipCPT.getColor((float)rawLen));
				}
				
				if (labelMagnitudes && i == result.vectorLocs.length-1) {
					// figure out parallel and perpendicular components
					double parallelComponent = VectorComponent.PARALLEL.getMagnitude(vector, lineAz);
					// Perpendicular component (subtract parallel component from the original vector)
					double perpendicularComponent = VectorComponent.PERPENDICULAR.getMagnitude(vector, lineAz);
					
					if (D) {
						System.out.println("Final Vector that ends at "+result.endLoc+": "+vector);
						System.out.println("\tLineAz: "+lineAz+"; dirV="+directionVector(lineAz));
						System.out.println("\tParallel: "+parallelComponent);
						System.out.println("\tPerpendicular: "+perpendicularComponent);
					}
					
					// arranged from top to bottom
					String[] labels = {
						"\u007C"+slipDF.format(rawLen)+"\u007C",
						angleDF.format(az)+"\u00B0",
						slipDF.format(parallelComponent)+" \u2225",
						slipDF.format(perpendicularComponent)+" \u22A5"
					};
					
					Preconditions.checkState(labels.length == annBufferDists.length);
					
					for (int l=0; l<labels.length; l++) {
						Location annLoc = locConstPlotDist(loc, centerLoc, lineAz, annBufferDists[l]);
						XYTextAnnotation ann = new XYTextAnnotation(labels[l], annLoc.lon, annLoc.lat);
						ann.setFont(numberAnnFont);
						ann.setTextAnchor(TextAnchor.CENTER);
						mapMaker.addAnnotation(ann);
					}
				}
			}
		}
		
//		if (!plotVectors && !integrals.isEmpty()) {
//			// draw line directions
//			double arrowLen = 0.2 * maxVectDist;
////			List<LineIntegralResult> arrowIntegrals = integrals.size() == 1 ?
////					List.of(integrals.get(0)) : List.of(integrals.get(0), integrals.get(integrals.size()-1));
//			List<LineIntegralResult> arrowIntegrals = integrals;
//			List<LocationList> arrows = new ArrayList<>(arrowIntegrals.size());
//			for (LineIntegralResult integral : arrowIntegrals) {
//				arrows.add(LocationList.of(integral.startLoc, integral.endLoc));
////				Location end = integral.endLoc;
////				double az = integral.azimuth;
////				Location arrowStart = locConstPlotDist(end, centerLoc, az - 135, arrowLen);
////				Location arrowEnd = locConstPlotDist(end, centerLoc, az + 135, arrowLen);
////				LocationList arrowhead = LocationList.of(arrowStart, end, arrowEnd);
////				lines.add(arrowhead);
////				lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//			}
////			mapMaker.plotArrows(arrows, arrowLen, Color.BL, 0)
//		}
		if (!lines.isEmpty()) {
//			System.out.println("Plotting "+lines.size()+" lines");
			mapMaker.plotLines(lines, lineChars);
		}
		if (!arrows.isEmpty()) {
//			System.out.println("Plotting "+arrows.size()+" arrows");
			double arrowLen;
			float arrowThickness;
			if (plotVectors) {
				arrowLen = 0.2*maxVectDist;
				arrowThickness = 4f;
			} else {
				arrowLen = 0.15*maxVectDist;
				arrowThickness = 2f;
			}
			mapMaker.setFillArrowheads(true);
			mapMaker.setArrowAngle(35);
			mapMaker.plotArrows(arrows, arrowLen, arrowColors, arrowThickness);
		}
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 5f);
		mapMaker.plotScatters(vectorLocs, Color.GRAY);
		
		return mapMaker;
	}
	
//	private static 
	
	public void plotIntegrals(File outputDir, String prefix, String title, boolean byLatitude,
			List<LineIntegralResult> integrals, VectorComponent... components) throws IOException {
		PlotSpec plot = buildIntegralPlot(title, byLatitude, integrals, components);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		if (components == null || components.length != 1) {
			// build our own ranges to leave room for the legend
			gp.drawGraphPanel(plot, false, false, null, getPlotYRange(plot));
		} else {
			gp.drawGraphPanel(plot);
		}
		
		PlotUtils.writePlots(outputDir, prefix, gp, 850, 800, true, true, false);
	}
	
	public static Range getPlotYRange(PlotSpec plot) {
		double maxY = 0d;
		double minY = Double.POSITIVE_INFINITY;
		for (DiscretizedFunc func : plot.getPlotFunctionsOnly()) {
			maxY = Math.max(maxY, func.getMaxY());
			minY = Math.min(minY, func.getMinY());
			if (func instanceof UncertainBoundedDiscretizedFunc) {
				maxY = Math.max(maxY, ((UncertainBoundedDiscretizedFunc)func).getUpperMaxY());
				minY = Math.min(minY, ((UncertainBoundedDiscretizedFunc)func).getLowerMinY());
			}
		}
		if (minY < 0 && maxY == 0d)
			// all below 0
			return new Range(Math.min(minY*1.2, minY - 3), 2);
		if (minY < 0)
			// below and above
			return new Range(Math.min(minY*1.2, minY - 3), Math.max(maxY*1.2, maxY + 3));
		// all above
		return new Range(0, Math.max(maxY*1.2, maxY + 3));
	}
	
	public DiscretizedFunc buildIntegralFunction(boolean byLatitude, List<LineIntegralResult> integrals, VectorComponent component) {
		DiscretizedFunc func = new ArbDiscrEmpiricalDistFunc();
		for (LineIntegralResult integral: integrals) {
			double x = byLatitude ? integral.endLoc.lat : integral.endLoc.lon;
			double y;
			if (integral.vectors.length == 0) {
				y = 0d;
			} else {
				Vector3D fullVect = integral.vectors[integral.vectors.length-1];
				y = component.getMagnitude(fullVect, integral.azimuth);
			}
			func.set(x, y);
		}
		return func;
	}
	
	public PlotSpec buildIntegralPlot(String title, boolean byLatitude, List<LineIntegralResult> integrals, VectorComponent... components) {
		if (components == null || components.length == 0)
			components = VectorComponent.values();
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int c=0; c<components.length; c++) {
			VectorComponent component = components[c];
			DiscretizedFunc func = buildIntegralFunction(byLatitude, integrals, component);
			
			func.setName(component.label);
			funcs.add(func);
			chars.add(new PlotCurveCharacterstics(components.length == 1 ? PlotLineType.SOLID : component.lineType, 3f, Color.BLACK));
		}
		
		String xAxisLabel = "Line End "+(byLatitude ? "Latitude" : "Longitude");
		String yAxisLabel = "Summed"+(components.length == 1 ? " "+components[0].label : "")+" Rate (mm/yr)";
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
//		spec.setLegendInset(components.length > 1);
		if (components.length > 1)
			spec.setLegendInset(RectangleAnchor.TOP_RIGHT, 0.975, 0.975, 0.8, false);
		
		return spec;
	}
	
	public void plotCombinedMapIntegrals(File outputDir, String prefix, String title, boolean byLatitude,
			List<LineIntegralResult> integrals, List<LineIntegralResult> mapIntegrals, boolean labelMagnitudes,
			boolean plotVectors, double cptMax, VectorComponent... components) throws IOException {
		GeographicMapMaker mapMaker = buildMapPlot(labelMagnitudes, plotVectors, cptMax, mapIntegrals);
		PlotSpec mapPlot = mapMaker.buildPlot(title);
		// move slip cpt legend to the top
		if (!byLatitude)
			for (Title subtitle : mapPlot.getSubtitles())
				subtitle.setPosition(RectangleEdge.TOP);
		PlotSpec integralPlot = buildIntegralPlot(title, byLatitude, integrals, components);
		
		Range lonRange = mapMaker.getXRange();
		Range latRange = mapMaker.getYRange();
		
		List<Range> xRanges = new ArrayList<>();
		List<Range> yRanges = new ArrayList<>();
		xRanges.add(lonRange);
		yRanges.add(latRange);
		
		List<PlotSpec> specs = new ArrayList<>();
		specs.add(mapPlot);
		
		if (byLatitude) {
			specs.add(getFlippedIntegralPlot(integralPlot));
			xRanges.add(getPlotYRange(integralPlot));
		} else {
			specs.add(integralPlot);
			yRanges.add(getPlotYRange(integralPlot));
		}

		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		PlotUtils.setSubPlotWeights(gp, getCombinedMapSubplotWeights(1, latRange, lonRange, byLatitude));
		
		PlotUtils.writePlots(outputDir, prefix, gp, byLatitude ? 1200 : 800, true, true, true, false);
	}
	
	public void plotCombinedMapMultiIntegrals(File outputDir, String prefix, String title, boolean byLatitude,
			List<List<LineIntegralResult>> integralsList, List<String> integralsLabels,
			List<LineIntegralResult> mapIntegrals, boolean sameSummedRange, double cptMax,
			VectorComponent... components) throws IOException {
		GeographicMapMaker mapMaker = buildMapPlot(false, false, cptMax, mapIntegrals);
		PlotSpec mapPlot = mapMaker.buildPlot(title);
		// move slip cpt legend to the top
		if (!byLatitude && cptMax > 0d)
			for (Title subtitle : mapPlot.getSubtitles())
				subtitle.setPosition(RectangleEdge.TOP);
		
		Range lonRange = mapMaker.getXRange();
		Range latRange = mapMaker.getYRange();
		
		List<Range> xRanges = new ArrayList<>();
		List<Range> yRanges = new ArrayList<>();
		xRanges.add(lonRange);
		yRanges.add(latRange);
		
		List<PlotSpec> specs = new ArrayList<>();
		specs.add(mapPlot);
		
		Preconditions.checkState(integralsList.size() == integralsLabels.size());
		
		for (int i=0; i<integralsList.size(); i++) {
			List<LineIntegralResult> integrals = integralsList.get(i);
			String label = integralsLabels.get(i);
			
			PlotSpec integralPlot = buildIntegralPlot(title, byLatitude, integrals, components);
			integralPlot.setYAxisLabel(label+" Rate (mm/yr)");
			
			Range range = getPlotYRange(integralPlot);
			
			if (byLatitude) {
				integralPlot = getFlippedIntegralPlot(integralPlot);
				xRanges.add(range);
			} else {
				yRanges.add(range);
			}
			
			if (i > 0)
				integralPlot.setLegendVisible(false);
			specs.add(integralPlot);
		}
		
		if (sameSummedRange) {
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;
			for (int i=1; i<specs.size(); i++) {
				Range range = byLatitude ? xRanges.get(i) : yRanges.get(i);
				min = Math.min(range.getLowerBound(), min);
				max = Math.max(range.getUpperBound(), max);
			}
			Range range = new Range(min, max);
			for (int i=1; i<specs.size(); i++) {
				if (byLatitude)
					xRanges.set(i, range);
				else
					yRanges.set(i, range);
			}
		}

		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		PlotUtils.setSubPlotWeights(gp, getCombinedMapSubplotWeights(integralsList.size(), latRange, lonRange, byLatitude));
		
		PlotUtils.writePlots(outputDir, prefix, gp, byLatitude ? 1200 : 800, true, true, true, false);
	}
	
	private static int[] getCombinedMapSubplotWeights(int numCharts, Range latRange, Range lonRange, boolean byLatitude) {
		double mapAspectRatio = PlotUtils.calcAspectRatio(lonRange, latRange, true);
		Preconditions.checkState(mapAspectRatio > 0d);
		double chartAspectRatio;
		if (numCharts >= 4)
			chartAspectRatio = 3d;
		else if (numCharts == 3)
			chartAspectRatio = 25d/10d;
		else if (numCharts == 2)
			chartAspectRatio = 18d/10d;
		else
			chartAspectRatio = 15d/10d;
		if (byLatitude)
			// chart(s) will be on the right, we want it taller than wide 
			chartAspectRatio = 1d / chartAspectRatio;
		double aspectScale = mapAspectRatio / chartAspectRatio;
		int mapWeight = 1000;
		int chartWeight = (int)((double)mapWeight * aspectScale + 0.5);
		int[] weights = new int[numCharts+1];
		weights[0] = mapWeight;
		for (int i=1; i<weights.length; i++)
			weights[i] = chartWeight;
		return weights;
	}
	
	private static PlotSpec getFlippedIntegralPlot(PlotSpec integralPlot) {
		List<XY_DataSet> flippedFuncs = new ArrayList<>();
		for (DiscretizedFunc func : integralPlot.getPlotFunctionsOnly()) {
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			for (Point2D pt : func)
				xy.set(pt.getY(), pt.getX());
			flippedFuncs.add(xy);
		}
		
		PlotSpec flipped = new PlotSpec(flippedFuncs, integralPlot.getChars(),
				integralPlot.getTitle(), integralPlot.getYAxisLabel(), integralPlot.getXAxisLabel());
		flipped.setLegendInset(integralPlot.isLegendInset());
		return flipped;
	}
	
	public void writeIntegralCSV(File outputFile, List<LineIntegralResult> integrals) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = List.of("Start Latitude", "Start Longitude", "End Latitude", "End Longitude",
				"Line Azimuth (degrees)", "Number of Faults Crossed", "Vector Azimuth (degrees)",
				"Summed Horizontal Rate (mm/yr)", "Line-Parallel Summed Rate (mm/yr)",
				"Line-Perpendicular Summed Rate (mm/yr)", "Vertical Summed Rate (mm/yr)");
		csv.addLine(header);
		
		DecimalFormat locDF = new DecimalFormat("0.#");
		DecimalFormat slipDF = new DecimalFormat("0.###");
		DecimalFormat azDF = new DecimalFormat("0");
		
		for (LineIntegralResult integral : integrals) {
			List<String> line = new ArrayList<>(header.size());
			line.add(locDF.format(integral.startLoc.lat));
			line.add(locDF.format(integral.startLoc.lon));
			line.add(locDF.format(integral.endLoc.lat));
			line.add(locDF.format(integral.endLoc.lon));
			line.add(azDF.format(integral.azimuth));
			int numCrossed = integral.vectors.length == 0 ? 0 : integral.vectors.length-2;
			line.add(numCrossed+"");
			if (numCrossed == 0) {
				while (line.size() < header.size())
					line.add("");
			} else {
				// last vector is summed across all
				Vector3D vect = integral.vectors[integral.vectors.length-1];
				line.add(azDF.format(vectorAzimuth(vect, integral.endLoc)));
				line.add(slipDF.format(VectorComponent.FULL_HORIZONTAL.getMagnitude(vect, integral.azimuth)));
				line.add(slipDF.format(VectorComponent.PARALLEL.getMagnitude(vect, integral.azimuth)));
				line.add(slipDF.format(VectorComponent.PERPENDICULAR.getMagnitude(vect, integral.azimuth)));
				line.add(slipDF.format(vect.getZ()));
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(outputFile);
	}
	
	public static Location locConstPlotDist(Location fromLoc, Location refCenterLoc, double azimuthDeg, double distance) {
		// calculate the length in the middle of the region so we don't get any distortion with large regions
		Location toLoc = LocationUtils.location(refCenterLoc, Math.toRadians(azimuthDeg), distance);
		// now move to the correct spot
		return new Location(fromLoc.lat - (refCenterLoc.lat - toLoc.lat), fromLoc.lon - (refCenterLoc.lon - toLoc.lon));
	}
	
	public static List<FaultSection> unionSubSectLists(List<List<? extends FaultSection>> sectsLists) {
		List<FaultSection> ret = new ArrayList<>();
		
		for (List<? extends FaultSection> sects : sectsLists) {
			for (FaultSection sect : sects) {
				sect = sect.clone();
				sect.setSectionId(ret.size());
				ret.add(sect);
			}
		}
		
		return ret;
	}
	
	private static Vector3D directionVector(double azimuth) {
		double azRad = Math.toRadians(azimuth);
		return new Vector3D(Math.sin(azRad), Math.cos(azRad), 0);
	}

	public static void main(String[] args) throws IOException {
//		double[] azimuths = { 0d, 45d, 90d, 135d, 180d, 225d, 270d, 315d, 360d };
//		double refAz = 15d;
//		for (double lineAz : azimuths) {
//			Vector3D lineDirection = directionVector(lineAz);
//			System.out.println("Azimuth: "+(float)lineAz+" deg");
//			System.out.println("\tDirection: "+lineDirection);
//			System.out.println("\tParallel (to "+(float)refAz+"): "+VectorComponent.PARALLEL.get(lineDirection, refAz));
//			System.out.println("\tPerpendicual (to "+(float)refAz+"): "+VectorComponent.PERPENDICULAR.get(lineDirection, refAz));
//		}
//		System.exit(0);
//		// fault 1
//		FaultTrace trace = new FaultTrace("test trace");
//		trace.add(new Location(0d, 0d));
//		trace.add(new Location(1d, 0d));
////		trace.add(new Location(1d, 1d));
////		double rake = 180d;
//		double rake = -90d;
////		double rake = -180d;
////		double dip = 90d;
//		double dip = 45d;
//		System.out.println("strike: "+trace.getAveStrike());
//
//		double upperDepth = 0d;
//		double lowerDepth = 10d;
//
//		FaultSectionPrefData sect = new FaultSectionPrefData();
//		sect.setAveDip(dip);
//		sect.setDipDirection((float)(90d + trace.getAveStrike()));
//		sect.setAveRake(rake);
//		sect.setFaultTrace(trace);
//		sect.setAveUpperDepth(upperDepth);
//		sect.setAveLowerDepth(lowerDepth);
//		
//		Vector3D vect = calcSlipVector(sect);
//		
//		System.out.println("Slip vector: "+vect+" [len="+(float)vectorLen(vect)+"]");
		
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/"
//				+ "2024_10_24-prvi25_subduction_branches/results_PRVI_SUB_FM_LARGE_branch_averaged.zip"));
//		
//		FaultSystemLineIntegralCalculator calc = new FaultSystemLineIntegralCalculator(sol.getRupSet());
		
		File outputDir = new File("/tmp/prvi_line_integrals");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<List<? extends FaultSection>> plotSects = new ArrayList<>();
		List<String> plotPrefixes = new ArrayList<>();
		List<String> plotTitles = new ArrayList<>();
		List<Double> plotMaxSlips = new ArrayList<>();
		List<List<List<? extends FaultSection>>> plotSectComponents = new ArrayList<>();
		
		List<? extends FaultSection> subLargeFullSects = PRVI25_SubductionDeformationModels.FULL.build(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE);
		List<? extends FaultSection> subSmallFullSects = PRVI25_SubductionDeformationModels.FULL.build(PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL);
		List<? extends FaultSection> subLargePartialSects = PRVI25_SubductionDeformationModels.PARTIAL.build(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE);
		List<? extends FaultSection> subSmallPartialSects = PRVI25_SubductionDeformationModels.PARTIAL.build(PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL);
//		List<? extends FaultSection> crustalSects = PRVI25_CrustalDeformationModels.GEOLOGIC.build(PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1);
		List<? extends FaultSection> crustalSects = PRVI25_CrustalDeformationModels.GEOLOGIC_DIST_AVG.build(PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1);

		plotSectComponents.add(List.of(subLargeFullSects, crustalSects));
		plotSects.add(unionSubSectLists(plotSectComponents.get(plotSectComponents.size()-1)));
		plotPrefixes.add("line_integral_sub_large_full_plus_crustal");
		plotTitles.add("Large Subduction, Full Rate + Crustal");
		plotMaxSlips.add(10d);
		
		plotSectComponents.add(List.of(subLargePartialSects, crustalSects));
		plotSects.add(unionSubSectLists(plotSectComponents.get(plotSectComponents.size()-1)));
		plotPrefixes.add("line_integral_sub_large_partial_plus_crustal");
		plotTitles.add("Large Subduction, Partial Rate + Crustal");
		plotMaxSlips.add(10d);
		
		plotSectComponents.add(List.of(subSmallFullSects, crustalSects));
		plotSects.add(unionSubSectLists(plotSectComponents.get(plotSectComponents.size()-1)));
		plotPrefixes.add("line_integral_sub_small_full_plus_crustal");
		plotTitles.add("Small Subduction, Full Rate + Crustal");
		plotMaxSlips.add(10d);
		
		plotSectComponents.add(List.of(subSmallPartialSects, crustalSects));
		plotSects.add(unionSubSectLists(plotSectComponents.get(plotSectComponents.size()-1)));
		plotPrefixes.add("line_integral_sub_small_partial_plus_crustal");
		plotTitles.add("Small Subduction, Partial Rate + Crustal");
		plotMaxSlips.add(10d);
		
		plotSectComponents.add(null);
		plotSects.add(subLargeFullSects);
		plotPrefixes.add("line_integral_sub_large_full_only");
		plotTitles.add("Large Subduction, Full Rate");
		plotMaxSlips.add(5d);
		
		plotSectComponents.add(null);
		plotSects.add(subLargePartialSects);
		plotPrefixes.add("line_integral_sub_large_partial_only");
		plotTitles.add("Large Subduction, Partial Rate");
		plotMaxSlips.add(5d);
		
		plotSectComponents.add(null);
		plotSects.add(subSmallFullSects);
		plotPrefixes.add("line_integral_sub_small_full_only");
		plotTitles.add("Small Subduction, Full Rate");
		plotMaxSlips.add(5d);
		
		plotSectComponents.add(null);
		plotSects.add(subSmallPartialSects);
		plotPrefixes.add("line_integral_sub_small_partial_only");
		plotTitles.add("Small Subduction, Partial Rate");
		plotMaxSlips.add(5d);
		
		plotSectComponents.add(null);
		plotSects.add(crustalSects);
		plotPrefixes.add("line_integral_crustal_only");
		plotTitles.add("Crustal");
		plotMaxSlips.add(10d);
		
		FaultSystemLineIntegralCalculator.D = false;
		
		double minLat = 16;
		double maxLat = 21;
		double minLon = -71d;
		double maxLon = -60.5d;
		double mapDelta = 1d;
		double chartDelta = 0.1;
		boolean byLatitude = false;

		List<FaultSystemLineIntegralCalculator> calcs = new ArrayList<>();
		List<List<LineIntegralResult>> mapResultsList = new ArrayList<>();
		List<List<LineIntegralResult>> denseResultsList = new ArrayList<>();
		
		boolean reducedSlipRates = true;
		
		for (int p=0; p<plotSects.size(); p++) {
			List<? extends FaultSection> sects = plotSects.get(p);
			String prefix = plotPrefixes.get(p);
			String title = plotTitles.get(p);
			double maxSlip = plotMaxSlips.get(p);
			
			System.out.println("Doing "+prefix+": "+title);
			
			FaultSystemLineIntegralCalculator calc = new FaultSystemLineIntegralCalculator(
					sects, reducedSlipRates);
			calcs.add(calc);
			
			boolean labelMags = true;
			System.out.println("Building map");
			List<LineIntegralResult> mapResults = new ArrayList<>();
			if (byLatitude) {
				for (double lat=minLat; (float)lat<(float)maxLat; lat += mapDelta)
					mapResults.add(calc.calcLineIntegral(new Location(lat, minLon), new Location(lat, maxLon)));
			} else {
				for (double lon=minLon; (float)lon<=(float)maxLon; lon += mapDelta)
					mapResults.add(calc.calcLineIntegral(new Location(minLat, lon), new Location(maxLat, lon)));
//				mapResults.add(calc.calcLineIntegral(new Location(minLat, -65.5), new Location(maxLat, -65.5)));
			}
			// prune empties on the ends
			while (mapResults.get(0).vectors.length == 0)
				mapResults.remove(0);
			while (mapResults.get(mapResults.size()-1).vectors.length == 0)
				mapResults.remove(mapResults.size()-1);
			mapResultsList.add(mapResults);
			System.out.println("\tHave "+mapResults.size()+" map integrals");
			
			calc.plotMap(outputDir, prefix+"_map", title, labelMags, true, maxSlip, mapResults);
			
			System.out.println("Building profiles");
			List<LineIntegralResult> denseResults = new ArrayList<>();
			if (byLatitude) {
				for (double lat=minLat; (float)lat<(float)maxLat; lat += chartDelta)
					denseResults.add(calc.calcLineIntegral(new Location(lat, minLon), new Location(lat, maxLon)));
			} else {
				for (double lon=minLon; (float)lon<=(float)maxLon; lon += chartDelta)
					denseResults.add(calc.calcLineIntegral(new Location(minLat, lon), new Location(maxLat, lon)));
			}
			denseResultsList.add(denseResults);
			System.out.println("\tHave "+denseResults.size()+" profile integrals");
			calc.plotIntegrals(outputDir, prefix, title, byLatitude, denseResults);
			
			calc.plotCombinedMapIntegrals(outputDir, prefix+"_comb", title, byLatitude, denseResults, mapResults, false, true, maxSlip);
			calc.plotCombinedMapIntegrals(outputDir, prefix+"_comb_no_vects", title, byLatitude, denseResults, mapResults, true, false, maxSlip);
			calc.writeIntegralCSV(new File(outputDir, prefix+".csv"), denseResults);
			
			List<List<? extends FaultSection>> components = plotSectComponents.get(p);
			if (components != null) {
				// now write out by component
				Preconditions.checkState(components.size() == 2);
				
				List<List<LineIntegralResult>> componentResultsList = new ArrayList<>();
				List<String> componentLabels = new ArrayList<>();
				
				componentResultsList.add(denseResults);
				componentLabels.add("Summed");
				
				componentLabels.add("Subduction");
				componentLabels.add("Crustal");
				
				for (List<? extends FaultSection> componentSects : components) {
					FaultSystemLineIntegralCalculator componentCalc = new FaultSystemLineIntegralCalculator(componentSects, reducedSlipRates) ;
					List<LineIntegralResult> componentResults = new ArrayList<>();
					for (LineIntegralResult integral : denseResults)
						componentResults.add(componentCalc.calcLineIntegral(integral.startLoc, integral.endLoc));
					componentResultsList.add(componentResults);
				}
				
				calc.plotCombinedMapMultiIntegrals(outputDir, prefix+"_by_component", title, byLatitude,
						componentResultsList, componentLabels, mapResults, true, maxSlip);
			}
		}
	}

}
