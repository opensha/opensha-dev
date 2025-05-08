package scratch.kevin.nshm23.devinSlipRateTests;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.faultSysSolution.util.minisections.MinisectionMappings;
import org.opensha.sha.earthquake.faultSysSolution.util.minisections.MinisectionSlipRecord;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

@Affects(FaultSystemRupSet.SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
@Affects(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
@Affects(FaultSystemSolution.RATES_FILE_NAME)
public enum DevinModDeformationModels implements RupSetDeformationModel {
	GEO_FROM_DEVIN("Geologic Rates Modified By Devin", "GeoModDevin", NSHM23_DeformationModels.GEOLOGIC) {

		@Override
		protected Map<Integer, SlipOverride> getSlipRateOverrides(List<? extends FaultSection> sects) throws IOException {
			CSVFile<String> csv = checkLoadCSV();
			Map<Integer, SlipOverride> ret = new HashMap<>();
			
			for (int row=1; row<csv.getNumRows(); row++) {
				int id = csv.getInt(row, 0);
				String name = csv.get(row, 1);
				double slipRate = csv.getDouble(row, 2);
				Preconditions.checkState(!ret.containsKey(id));
				ret.put(id, new SlipOverride(id, name, slipRate));
			}
			
			return ret;
		}
		
	},
	GEO_AVG_FROM_DEVIN("Geologic Rates Modified By Devin, Averaged Per Parent", "GeoAvgModDevin", NSHM23_DeformationModels.GEOLOGIC) {

		@Override
		protected Map<Integer, SlipOverride> getSlipRateOverrides(List<? extends FaultSection> sects) throws IOException {
			Map<Integer, SlipOverride> origByID = GEO_FROM_DEVIN.getSlipRateOverrides(null);
			
			// aggregate by parent
			Map<Integer, List<SlipOverride>> byParent = new HashMap<>();
			
			for (SlipOverride slip : origByID.values()) {
				FaultSection sect = sects.get(slip.id);
				Preconditions.checkState(sect.getName().equals(slip.name));
				int parentID = sect.getParentSectionId();
				List<SlipOverride> parentSlips = byParent.get(parentID);
				if (parentSlips == null) {
					parentSlips = new ArrayList<>();
					byParent.put(parentID, parentSlips);
				}
				parentSlips.add(slip);
			}
			
			// average across each parent
			Map<Integer, SlipOverride> ret = new HashMap<>();
			
			for (int parentID : byParent.keySet()) {
				List<SlipOverride> slips = byParent.get(parentID);
				double avgSlip = slips.stream().mapToDouble(S->S.slipRate).average().getAsDouble();
				for (SlipOverride slip : slips)
					ret.put(slip.id, new SlipOverride(slip.id, slip.name, avgSlip));
			}
			
			return ret;
		}
		
	};
	
	static class SlipOverride {
		final int id;
		final String name;
		final double slipRate;
		private SlipOverride(int id, String name, double slipRate) {
			super();
			this.id = id;
			this.name = name;
			this.slipRate = slipRate;
		}
	}

	private String name;
	private String shortName;
	private NSHM23_DeformationModels refDM;

	private DevinModDeformationModels(String name, String shortName, NSHM23_DeformationModels refDM) {
		this.name = name;
		this.shortName = shortName;
		this.refDM = refDM;
	}
	
	private static String CSV_PATH = "/scratch/kevin/nshm23/devinSlipRateTests/devin_slip_modifications.csv";
	private static CSVFile<String> cachedCSV;
	
	abstract Map<Integer, SlipOverride> getSlipRateOverrides(List<? extends FaultSection> sects) throws IOException;
	
	
	private static CSVFile<String> checkLoadCSV() throws IOException {
		if (cachedCSV != null)
			return cachedCSV;
		synchronized (DevinModDeformationModels.class) {
			if (cachedCSV == null) {
				InputStream is = DevinModDeformationModels.class.getResourceAsStream(CSV_PATH);
				Preconditions.checkNotNull(is, "Null input stream for %s", CSV_PATH);
				cachedCSV = CSVFile.readStream(is, true);
			}
		}
		return cachedCSV;
	}

	@Override
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return 1d;
	}

	@Override
	public String getFilePrefix() {
		return name();
	}

	@Override
	public String getShortName() {
		return shortName;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public boolean isApplicableTo(RupSetFaultModel faultModel) {
		return faultModel instanceof NSHM23_FaultModels;
	}
		
//	@Override
//	public List<? extends FaultSection> build(RupSetFaultModel faultModel) throws IOException {
//		return build(faultModel, MIN_SUB_SECTS_PER_FAULT_DEFAULT, DOWN_DIP_FRACT_DEFAULT, MAX_LEN_DEFAULT);
//	}
//	
//	@Override
//	public List<? extends FaultSection> build(RupSetFaultModel faultModel, int minPerFault, double ddwFract,
//			double fixedLen) throws IOException {
//		Map<Integer, List<MinisectionSlipRecord>> minis = refDM.getMinisections(faultModel);
//		
//		// fetch full sections
//		List<? extends FaultSection> fullSects = faultModel.getFaultSections();
//
//		// no subsections passed in, build them
//		List<? extends FaultSection> subSects = SubSectionBuilder.buildSubSects(
//				fullSects, minPerFault, ddwFract, fixedLen);
//
//		return buildDeformationModel(faultModel, minis, fullSects, subSects);
//	}
//
//	@Override
//	public List<? extends FaultSection> buildForSubsects(
//			RupSetFaultModel faultModel, List<? extends FaultSection> subSects) throws IOException {
//		Map<Integer, List<MinisectionSlipRecord>> minis = refDM.getMinisections(faultModel);
//		
//		// fetch full sections
//		List<? extends FaultSection> fullSects = faultModel.getFaultSections();
//		
//		return buildDeformationModel(faultModel, minis, fullSects, subSects);
//	}
	
	@Override
	public List<? extends FaultSection> apply(RupSetFaultModel faultModel,
			LogicTreeBranch<? extends LogicTreeNode> branch, List<? extends FaultSection> fullSects,
			List<? extends FaultSection> subSects) throws IOException {
		Map<Integer, List<MinisectionSlipRecord>> minis = refDM.getMinisections(faultModel);
		
		// minisection mappings
		MinisectionMappings mappings = new MinisectionMappings(fullSects, subSects);
		
		List<? extends FaultSection> origSects = refDM.apply(faultModel, branch, fullSects, subSects);
		Preconditions.checkState(origSects.size() == subSects.size());
		
		System.out.println("Building modified deformation model with "+subSects.size()+" sub-sections: "+name);
		Map<Integer, SlipOverride> modSlipRates = getSlipRateOverrides(subSects);
		
		for (int i=0; i<subSects.size(); i++) {
			FaultSection subSect = subSects.get(i);
			FaultSection origSect = origSects.get(i);
			Preconditions.checkState(subSect.getSectionName().equals(origSect.getSectionName()));
			double origSlipRate = origSect.getOrigAveSlipRate();
			if (modSlipRates.containsKey(i)) {
//			if ((float)modSlipRates[i] != (float)origSlipRate) {
				SlipOverride modSlipRate = modSlipRates.get(i);
				Preconditions.checkState(modSlipRate.name.equals(subSect.getName()),
						"Name mismatch for %s: '%s' != '%s'", i, subSect.getName(), modSlipRate.name);;
				System.out.println(i+". "+subSect.getSectionName()+":\t"+(float)origSlipRate+"\t->\t"+(float)modSlipRate.slipRate);
				subSect.setAveSlipRate(modSlipRate.slipRate);
			} else {
//				System.out.println("IDENTICAL "+i+". "+subSect.getSectionName()+":\t"+(float)origSlipRate);
				subSect.setAveSlipRate(origSlipRate);
			}
		}
		
		// apply std dev defaults
		return refDM.applyCreepModel(mappings,
				refDM.applyStdDevDefaults(faultModel, subSects));
	}


	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		for (DevinModDeformationModels dm : DevinModDeformationModels.values()) {
			dm.build(fm, null);
		}
	}

}
