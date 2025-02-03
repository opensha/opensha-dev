package scratch.kevin.prvi25;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

public class SubductionConvertMuertosToCrustal {

	public static void main(String[] args) throws IOException {
		ModuleContainer.VERBOSE_DEFAULT = false;
		File inputFile, outputFile;
		if (args.length == 0) {
			inputFile = COMBINED_SOL;
			outputFile = new File(COMBINED_DIR, "combined_branch_averaged_solution_mue_as_crustal.zip");
		} else {
			inputFile = new File(args[0]);
			outputFile = new File(args[1]);
		}
		
		FaultSystemSolution inputSol = FaultSystemSolution.load(inputFile);
		
		FaultSystemRupSet rupSet = inputSol.getRupSet();
		RupSetTectonicRegimes regimes = rupSet.requireModule(RupSetTectonicRegimes.class);
		
		TectonicRegionType[] modRegimes = new TectonicRegionType[rupSet.getNumRuptures()];
		for (int r=0; r<modRegimes.length; r++)
			modRegimes[r] = regimes.get(r);
		
		// now replace muertos fault based
		int mueParent = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Muertos");
		int numSwitched = 0;
		for (int rupIndex : rupSet.getRupturesForParentSection(mueParent)) {
			Preconditions.checkState(modRegimes[rupIndex] == TectonicRegionType.SUBDUCTION_INTERFACE);
			modRegimes[rupIndex] = TectonicRegionType.ACTIVE_SHALLOW;
			numSwitched++;
		}
		System.out.println("Switched "+numSwitched+" supra-seis Muertos ruptures from interface to active");
		
		rupSet.addModule(new RupSetTectonicRegimes(rupSet, modRegimes));
		
		// now replace muertos gridded
		GridSourceList origGridSources = inputSol.requireModule(GridSourceList.class);
		int numLocs = origGridSources.getNumLocations();
		
		EnumMap<TectonicRegionType, List<List<GriddedRupture>>> trtRuptureLists = new EnumMap<>(TectonicRegionType.class);
		for (TectonicRegionType trt : origGridSources.getTectonicRegionTypes()) {
			List<List<GriddedRupture>> list = new ArrayList<>(numLocs);
			for (int l=0; l<numLocs; l++)
				list.add(new ArrayList<>());
			trtRuptureLists.put(trt, list);
		}
		
		List<List<GriddedRupture>> outShallowRupLists = trtRuptureLists.get(TectonicRegionType.ACTIVE_SHALLOW);
		
		HashSet<FaultSection> allInterfaceAssocSects = new HashSet<>();
		
		int numGridSwitched = 0;
		for (int l=0; l<origGridSources.getNumLocations(); l++) {
			for (TectonicRegionType trt : TectonicRegionType.values()) {
				ImmutableList<GriddedRupture> rups = origGridSources.getRuptures(trt, l);
				if (!rups.isEmpty()) {
					List<GriddedRupture> outRups = trtRuptureLists.get(trt).get(l);
					for (GriddedRupture rup : rups) {
						if (rup.properties.tectonicRegionType == TectonicRegionType.SUBDUCTION_INTERFACE) {
							// can be one or more (multiple instances of same subsection)
							Preconditions.checkState(rup.associatedSections != null && rup.associatedSections.length >= 1,
									"Expected 1+ association for interface, have %s; rup: %s",
									rup.associatedSections == null ? 0 : rup.associatedSections.length, rup);
							boolean muertos = false;
							for (int i=0; i<rup.associatedSections.length; i++) {
								FaultSection sect = rupSet.getFaultSectionData(rup.associatedSections[0]);
								boolean myMuertos = sect.getParentSectionId() == mueParent;
								if (!allInterfaceAssocSects.contains(sect)) {
									System.out.println("Encountered interface associated subsect: "+sect.getName()+" (muertos="+myMuertos+")");
									allInterfaceAssocSects.add(sect);
								}
								if (i == 0)
									muertos = myMuertos;
								else
									Preconditions.checkState(muertos == myMuertos);
							}
							
							if (muertos) {
								// convert to crustal
								GriddedRuptureProperties props = rup.properties;
								GriddedRuptureProperties modProps = new GriddedRuptureProperties(
										props.magnitude, props.rake, props.dip, props.strike, props.strikeRange,
										props.upperDepth, props.lowerDepth, props.length, props.hypocentralDepth,
										props.hypocentralDAS, TectonicRegionType.ACTIVE_SHALLOW);
								outShallowRupLists.get(l).add(new GriddedRupture(l, rup.location, modProps,
										rup.rate, rup.associatedSections, rup.associatedSectionFracts));
								numGridSwitched++;
							} else {
								outRups.add(rup);
							}
						} else {
							outRups.add(rup);
						}
					}
				}
			}
		}
		System.out.println("Switched "+numGridSwitched+" gridded Muertos ruptures from interface to active");
		
		inputSol.addModule(new GridSourceList.Precomputed(origGridSources.getGriddedRegion(), trtRuptureLists));
		
		inputSol.write(outputFile);
	}

}