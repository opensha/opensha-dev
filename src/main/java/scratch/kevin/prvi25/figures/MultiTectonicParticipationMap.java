package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.jfree.chart.title.PaintScaleLegend;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class MultiTectonicParticipationMap {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "combined_partic_rate_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdirs());
		
		FaultSystemSolution crustalSol = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED);
		FaultSystemSolution subSmallSol = FaultSystemSolution.load(SUBDUCTION_SOL_SMALL);
		FaultSystemSolution subLargeSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		FaultSystemSolution combSol = FaultSystemSolution.load(COMBINED_SOL);
		
		boolean includeInterfaceGridded = false;
		
		double cptMin = -6d;
		double cptMax = -2d;
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(cptMin, cptMax);
//		CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(cptMin, cptMax);
//		CPT cpt = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(cptMin, cptMax);
//		CPT cpt = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(cptMin, cptMax);
//		CPT cpt = GMT_CPT_Files.SEQUENTIAL_OSLO_UNIFORM.instance().rescale(cptMin, cptMax);
//		CPT cpt = new CPT(cptMin, cptMax, new Color(0, 0, 255, 40), new Color(255, 0, 0, 40));
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setNanColor(Color.GRAY);
		cpt.setPreferredTickInterval(0.5);
		
		GridSourceList subSmallGridded = null;
		GridSourceList subLargeGridded = null;
		if (includeInterfaceGridded) {
			subSmallGridded = subSmallSol.requireModule(GridSourceList.class);
			subLargeGridded = subLargeSol.requireModule(GridSourceList.class);
		}
		
		List<? extends FaultSection> crustalSects = crustalSol.getRupSet().getFaultSectionDataList();
		List<? extends FaultSection> subSmallSects = subSmallSol.getRupSet().getFaultSectionDataList();
		List<? extends FaultSection> subLargeSects = subLargeSol.getRupSet().getFaultSectionDataList();
		System.out.println("SubSmall has "+subSmallSects.size()+" sects");
		System.out.println("SubLarge has "+subLargeSects.size()+" sects");
		
		boolean allMatch = subSmallSects.size() == subLargeSects.size();
		for (int i=0; i<Integer.max(subSmallSects.size(), subLargeSects.size()); i++) {
			String name1;
			if (i >= subSmallSects.size()) {
				name1 = "MISSING";
			} else {
				FaultSection sect = subSmallSects.get(i);
				name1 = sect.getSectionName()+" ("+sect.getParentSectionId()+")";
			}
			String name2;
			if (i >= subLargeSects.size()) {
				name2 = "MISSING";
			} else {
				FaultSection sect = subLargeSects.get(i);
				name2 = sect.getSectionName()+" ("+sect.getParentSectionId()+")";
			}
			boolean match = name1.equals(name2);
			System.out.println(i+". "+name1+"\t"+name2+(match ? "" : "\tMISMATCH"));
			allMatch &= match;
		}
		Preconditions.checkState(allMatch, "Sub FM sects don't match");

		double[] minMags = { 0d, 6d, 6.5, 7d, 7.5d, 8d };
		double[] gridMinMags = { 5d, 6d, 6.5, 7d, 7.5d, 8d };
		
		List<FaultSection> combSects = new ArrayList<>();
		List<List<Double>> combSectRates = new ArrayList<>();
		List<List<Double>> combSectSortables = new ArrayList<>();
		for (int m=0; m<minMags.length; m++) {
			combSectRates.add(new ArrayList<>());
			combSectSortables.add(new ArrayList<>());
		}
		
		// first crustal
		for (FaultSection sect : crustalSects) {
			for (int m=0; m<minMags.length; m++) {
				double rate = crustalSol.calcParticRateForSect(sect.getSectionId(), minMags[m], Double.POSITIVE_INFINITY);
				combSectRates.get(m).add(rate);
				combSectSortables.get(m).add(rate);
			}
			sect = sect.clone();
			sect.setSectionId(combSects.size());
			combSects.add(sect);
		}
		
		int muertosParent = FaultSectionUtils.findParentSectionID(subSmallSects, "Muertos");
		int hispaniolaParent = FaultSectionUtils.findParentSectionID(subSmallSects, "Hispaniola");
		double smallWeight = PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL.getNodeWeight(null);
		double largeWeight = PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE.getNodeWeight(null);
		if (smallWeight + largeWeight != 1d) {
			double sum = smallWeight + largeWeight;
			smallWeight /= sum;
			largeWeight /= sum;
		}
		
		List<double[]> smallGriddedRates = null;
		List<double[]> largeGriddedRates = null;
		if (includeInterfaceGridded) {
			for (boolean large : new boolean[] {false,true}) {
				List<double[]> griddedRates = new ArrayList<>(subSmallSects.size());
				for (int s=0; s<subSmallSects.size(); s++)
					griddedRates.add(new double[minMags.length]);
				GridSourceList grid = large ? subLargeGridded : subSmallGridded;
				for (int l=0; l<grid.getNumLocations(); l++) {
					for (GriddedRupture rup : grid.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, l)) {
						Preconditions.checkState(rup.associatedSections != null);
						for (int i=0; i<rup.associatedSections.length; i++) {
							double[] rates = griddedRates.get(rup.associatedSections[i]);
							for (int m=0; m<minMags.length; m++)
								if (rup.properties.magnitude >= minMags[m])
									rates[m] += rup.associatedSectionFracts[i] * rup.rate;
						}
					}
				}
				if (large)
					largeGriddedRates = griddedRates;
				else
					smallGriddedRates = griddedRates;
			}
		}
		
		// now subduction combined
		for (int s=0; s<subSmallSects.size(); s++) {
			FaultSection smallSect = subSmallSects.get(s);
			FaultSection largeSect = subLargeSects.get(s);
			
			double[] smallRates = null;
			double[] largeRates = null;
			
			for (boolean large : new boolean[] {false,true}) {
				double[] rates = new double[minMags.length];
				FaultSystemSolution sol = large ? subLargeSol : subSmallSol;
				for (int m=0; m<minMags.length; m++)
					rates[m] = sol.calcParticRateForSect(s, minMags[m], Double.POSITIVE_INFINITY);
				if (includeInterfaceGridded) {
					double[] gridRates = large ? largeGriddedRates.get(s) : smallGriddedRates.get(s);
					for (int m=0; m<minMags.length; m++)
						if (minMags[m] > 0d)
							rates[m] += gridRates[m];
				}
				
				if (large)
					largeRates = rates;
				else
					smallRates = rates;
			}
			
			if (smallSect.getParentSectionId() == muertosParent || smallSect.getParentSectionId() == hispaniolaParent) {
				// simple, just weight-average them and represent once
				for (int m=0; m<minMags.length; m++) {
					combSectRates.get(m).add(smallRates[m]*smallWeight + largeRates[m]*largeWeight);
					combSectSortables.get(m).add(-1d);
				}
				
				smallSect = smallSect.clone();
				smallSect.setSectionId(combSects.size());
				combSects.add(smallSect);
			} else {
				// put large on bottom with half weight, but small on top with full average
				for (int m=0; m<minMags.length; m++) {
					combSectRates.get(m).add(smallRates[m]*smallWeight + largeRates[m]*largeWeight);
					combSectSortables.get(m).add(-1d);
					
					combSectRates.get(m).add(largeRates[m]*largeWeight);
					combSectSortables.get(m).add(-2d);
				}
				
				smallSect = smallSect.clone();
				smallSect.setSectionId(combSects.size());
				smallSect.setSectionName(smallSect.getSectionName()+" (SMALL)");
				combSects.add(smallSect);

				largeSect = largeSect.clone();
				largeSect.setSectionId(combSects.size());
				largeSect.setSectionName(largeSect.getSectionName()+" (LARGE)");
				combSects.add(largeSect);
			}
		}
		
		// set zeros to NaNs
		for (List<Double> rates : combSectRates)
			for (int r=0; r<rates.size(); r++)
				if (rates.get(r) == 0d)
					rates.set(r, Double.NaN);
		
		// convert to log10
		for (List<Double> rates : combSectRates)
			for (int r=0; r<rates.size(); r++)
				rates.set(r, Math.log10(rates.get(r)));
		
		Region mapReg = PRVI25_RegionLoader.loadPRVI_ModelBroad();
		GeographicMapMaker mapMaker = new GeographicMapMaker(mapReg);
		
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setWritePDFs(true);
		mapMaker.setFillSurfaces(true);
		mapMaker.setSectOutlineChar(null);
		mapMaker.setReverseSort(false);
		mapMaker.setAbsoluteSort(false);
		mapMaker.setSectNaNChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, new Color(160, 160, 160, 140)));
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.DARK_GRAY));
		mapMaker.setScalarThickness(4f);
		mapMaker.setFaultSections(combSects);
		
		// write generic cpt only first
		PlotUtils.writeScaleLegendOnly(outputDir, "participation_cpt",
				GeographicMapMaker.buildCPTLegend(cpt, "Log10 Participation Rate (/yr)"),
				mapMaker.getDefaultPlotWidth(), true, true);
		PlotUtils.writeScaleLegendOnly(outputDir, "nucleation_cpt",
				GeographicMapMaker.buildCPTLegend(cpt, "Log10 Nucleation Rate (/yr)"),
				mapMaker.getDefaultPlotWidth(), true, true);
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		for (int i=0; i<combSects.size(); i++) {
			System.out.println(i+". "+combSects.get(i).getName()+":\trate="+combSectRates.get(0).get(i).floatValue()+":\tsort="+combSectSortables.get(0).get(i).floatValue());
		}
		
		for (int m=0; m<minMags.length; m++) {
			String label = "Log10 ";
			String prefix = "multi_tect_partic_";
			String noLegendTitle;
			if (minMags[m] > 0d) {
				label += "M>"+oDF.format(minMags[m]);
				prefix += "m"+(float)minMags[m];
				noLegendTitle = "Fault, M>"+oDF.format(minMags[m]);
			} else {
				label += "Supra-Seismogenic";
				prefix += "supra_seis";
				noLegendTitle = "Fault, Supra-Seismogenic";
			}
			if (includeInterfaceGridded)
				prefix += "_incl_interface_gridded";
			label += " Fault Participation Rate (/yr)";
			
			mapMaker.plotSectScalars(combSectRates.get(m), combSectSortables.get(m), cpt, label);
			mapMaker.plot(outputDir, prefix, " ");
			
			// now again without a CPT legend
			mapMaker.plotSectScalars(combSectRates.get(m), combSectSortables.get(m), cpt, null);
			mapMaker.plot(outputDir, prefix+"_no_cpt", noLegendTitle);
		}
		
		// now gridded
		TectonicRegionType[] trts = {
				null,
				TectonicRegionType.ACTIVE_SHALLOW,
				TectonicRegionType.SUBDUCTION_INTERFACE,
				TectonicRegionType.SUBDUCTION_SLAB
		};
		
		mapMaker.clearSectScalars();
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
		GridSourceList gridProv = combSol.requireModule(GridSourceList.class);
		GriddedRegion gridReg = gridProv.getGriddedRegion();
		cpt.setNanColor(Color.WHITE);
		double minGridMag = Doubles.min(gridMinMags);
		int numRandPartic = 1000;
		BitSet nodeBits = null;
		Table<Double, GriddedRuptureProperties, CachedGriddedParticResult> particsCache = HashBasedTable.create();
		Random rand = new Random(123456789l);
		for (boolean psuedoPartic : new boolean[] {false, true}) {
			if (psuedoPartic) {
				// precache
				ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
				Table<Double, GriddedRuptureProperties, Future<CachedGriddedParticResult>> particsCacheFutures = HashBasedTable.create();
				System.out.println("Precaching gridded participation maps");
				for (int l=0; l<gridProv.getNumLocations(); l++) {
					for (GriddedRupture rup : gridProv.getRuptures(null, l)) {
						if (rup.properties.magnitude < minGridMag || rup.properties.length == 0d || Double.isFinite(rup.properties.strike))
							continue;
						Location gridLoc = gridProv.getLocation(l);
						Future<CachedGriddedParticResult> future = particsCacheFutures.get(gridLoc.lat, rup.properties);
						if (future == null) {
							long seed = rand.nextLong();
							future = exec.submit(new Callable<CachedGriddedParticResult>() {

								@Override
								public CachedGriddedParticResult call() throws Exception {
									return new CachedGriddedParticResult(rup.properties, gridReg, gridLoc.lat, numRandPartic, new Random(seed));
								}
							});
							particsCacheFutures.put(gridLoc.lat, rup.properties, future);
						}
					}
				}
				System.out.println("Waiting on "+particsCacheFutures.size()+" futures");
				for (Cell<Double, GriddedRuptureProperties, Future<CachedGriddedParticResult>> cell : particsCacheFutures.cellSet()) {
					try {
						CachedGriddedParticResult result = cell.getValue().get();
						particsCache.put(cell.getRowKey(), cell.getColumnKey(), result);
					} catch (InterruptedException | ExecutionException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
				}
				
				exec.shutdown();
			}
			for (TectonicRegionType trt : trts) {
				GriddedGeoDataSet[] xyzs = new GriddedGeoDataSet[gridMinMags.length];
				for (int m=0; m<gridMinMags.length; m++)
					xyzs[m] = new GriddedGeoDataSet(gridReg);
				System.out.println("Doing gridded for partic="+psuedoPartic+", trt="+trt);
				for (int l=0; l<gridProv.getNumLocations(); l++) {
					for (GriddedRupture rup : gridProv.getRuptures(trt, l)) {
						if (rup.properties.magnitude < minGridMag)
							continue;
						if (psuedoPartic && rup.properties.length > 0d) {
							Location gridLoc = gridProv.getLocation(l);
							if (Double.isFinite(rup.properties.strike)) {
								// one off
								PointSurfaceBuilder builder = GridSourceList.surfBuilderForRup(rup);
								EvenlyGriddedSurface surf = rup.properties.dip == 90d ? builder.buildLineSurface() : builder.buildGriddedSurface();
								if (nodeBits == null)
									nodeBits = new BitSet(xyzs[0].size());
								else
									nodeBits.clear();
								for (Location loc : surf) {
									int index = gridReg.indexForLocation(loc);
									if (index >= 0)
										nodeBits.set(index);
								}
								for (int i=0; i<xyzs[0].size(); i++)
									if (nodeBits.get(i))
										for (int m=0; m<gridMinMags.length; m++)
											if (rup.properties.magnitude >= gridMinMags[m])
												xyzs[m].set(i, xyzs[m].get(i) + rup.rate);
							} else {
								// cached
								CachedGriddedParticResult result = particsCache.get(gridLoc.lat, rup.properties);
								if (result == null) {
									result = new CachedGriddedParticResult(rup.properties, gridReg, gridLoc.lat, numRandPartic, new Random(rand.nextLong()));
									particsCache.put(gridLoc.lat, rup.properties, result);
								}
								for (int xIndex=0; xIndex<result.xyz.getNumX(); xIndex++) {
									double lon = gridLoc.lon + result.xyz.getX(xIndex);
									for (int yIndex=0; yIndex<result.xyz.getNumY(); yIndex++) {
										double lat = gridLoc.lat + result.xyz.getY(yIndex);
										int index = gridReg.indexForLocation(new Location(lat, lon));
										if (index >= 0)
											for (int m=0; m<gridMinMags.length; m++)
												if (rup.properties.magnitude >= gridMinMags[m])
													xyzs[m].set(index, xyzs[m].get(index) + rup.rate*result.xyz.get(xIndex, yIndex));
									}
								}
							}
						} else {
							for (int m=0; m<gridMinMags.length; m++)
								if (rup.properties.magnitude >= gridMinMags[m])
									xyzs[m].set(l, xyzs[m].get(l)+rup.rate);
						}
					}
				}

				for (int m=0; m<gridMinMags.length; m++) {
					for (int l=0; l<xyzs[m].size(); l++)
						if (xyzs[m].get(l) == 0d)
							xyzs[m].set(l, Double.NaN);

					xyzs[m].log10();

					String noLegendTitle = "Gridded";
					String label = "Log10 ";
					String prefix = "gridded_";
					label += "M>"+oDF.format(gridMinMags[m]);
					prefix += "m"+(float)gridMinMags[m];
					if (trt != null) {
						switch (trt) {
						case ACTIVE_SHALLOW:
							label += " Crustal";
							noLegendTitle += ", Crustal";
							prefix += "_crustal";
							break;
						case SUBDUCTION_INTERFACE:
							label += " Interface";
							noLegendTitle += ", Interface";
							prefix += "_interface";
							break;
						case SUBDUCTION_SLAB:
							label += " Slab";
							noLegendTitle += ", Slab";
							prefix += "_slab";
							break;

						default:
							throw new IllegalStateException();
						}
					}
					noLegendTitle += ", M>"+oDF.format(gridMinMags[m]);
					if (psuedoPartic) {
						prefix += "_partic";
						label += " Gridded Participation Rate (/yr)";
					} else {
						label += " Gridded Nucleation Rate (/yr)";
					}

					mapMaker.plotXYZData(xyzs[m], cpt, label);
					mapMaker.plot(outputDir, prefix, " ");
					
					if (trt == null) {
						// now again without a CPT legend
						mapMaker.plotXYZData(xyzs[m], cpt, null);
						mapMaker.plot(outputDir, prefix+"_no_cpt", noLegendTitle);
					}
				}
			}
		}
	}
	
	private static class CachedGriddedParticResult {
		
		private EvenlyDiscrXYZ_DataSet xyz;
		private GriddedRegion surroundingGrid;

		public CachedGriddedParticResult(GriddedRuptureProperties props, GriddedRegion gridReg, double lat, int numSurfs, Random rand) {
			Preconditions.checkState(props.length > 0d);
			
			double regWidth = props.length;
			if (props.dip != 90d) {
				// sin(dip) = horz / vert
				// horz = vert * sin(dip)
				double ddwProj = Math.sin(Math.toRadians(props.dip)) * (props.lowerDepth - props.upperDepth);
				regWidth += ddwProj; // rotations get funky, buffer by some extra
			}
			regWidth += 50;
			
			Location refLoc = new Location(lat, 0.5*(gridReg.getMinLon() + gridReg.getMaxLon()));
			Location origRefLoc = refLoc;
			refLoc = gridReg.getLocation(gridReg.indexForLocation(refLoc)); // snap to grid
			Preconditions.checkNotNull(refLoc, "Couldn't snap middle loc to grid? orig=[%s]", origRefLoc);
			
			Preconditions.checkState(gridReg.getLatSpacing() == gridReg.getLonSpacing());
			double spacing = gridReg.getSpacing();
			
			Location bottom = LocationUtils.location(refLoc, Math.PI, 0.5*regWidth);
			Location bottomLeft = LocationUtils.location(bottom, 1.5*Math.PI, 0.5*regWidth);
			Location top = LocationUtils.location(refLoc, 0d, 0.5*regWidth);
			Location topRight = LocationUtils.location(top, 0.5*Math.PI, 0.5*regWidth);
			Region surroundingRegion = new Region(bottomLeft, topRight);
			
			surroundingGrid = new GriddedRegion(surroundingRegion, spacing, gridReg.getLocation(0));
			xyz = new EvenlyDiscrXYZ_DataSet(
					surroundingGrid.getNumLonNodes(), surroundingGrid.getNumLatNodes(),
					surroundingGrid.getMinGridLon()-refLoc.lon, surroundingGrid.getMinGridLat()-refLoc.lat, spacing);
			Preconditions.checkState(xyz.size() == surroundingGrid.getNodeCount());
			
			GriddedRupture fakeRup = new GriddedRupture(0, refLoc, props, 1d);
			
			PointSurfaceBuilder builder = GridSourceList.surfBuilderForRup(fakeRup);
			
			builder.random(rand);
			
			EvenlyGriddedSurface[] surfs;
			if (props.dip == 90d)
				surfs = builder.buildRandLineSurfaces(numSurfs);
			else
				surfs = builder.buildRandGriddedSurfaces(numSurfs);
			
			BitSet mappings = new BitSet(xyz.size());
			double rateEach = 1d/surfs.length;
			for (int s=0; s<surfs.length; s++) {
				if (s > 0)
					mappings.clear();
				for (Location loc : surfs[s]) {
					int mappedIndex = surroundingGrid.indexForLocation(loc);
					Preconditions.checkState(mappedIndex >= 0,
							"Buffer wasn't big enough? Surf location (%s) didn't map to index. Len=%s, regWidth=%s"
							+ "\n\trefLoc=[%s], botLeft=[%s], topRight=[%s]"
							+ "\n\tgridMin=[%s, %s], gridMax=[%s, %s]"
							+ "\n\tprops=%s",
							loc, props.length, regWidth, refLoc, bottomLeft, topRight,
							surroundingGrid.getMinGridLat(), surroundingGrid.getMinGridLon(),
							surroundingGrid.getMaxGridLat(), surroundingGrid.getMaxGridLon(),
							props);
					mappings.set(mappedIndex);
				}
				
				for (int i=0; i<xyz.size(); i++) {
					if (mappings.get(i)) {
						// don't assume order is same in geo and xyz
						Location loc = surroundingGrid.getLocation(i);
						int xIndex = xyz.getXIndex(loc.lon - refLoc.lon);
						int yIndex = xyz.getYIndex(loc.lat - refLoc.lat);
						xyz.set(xIndex, yIndex, xyz.get(xIndex, yIndex)+rateEach);
					}
				}
			}
		}
	}

}
