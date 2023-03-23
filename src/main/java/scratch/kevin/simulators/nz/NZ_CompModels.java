package scratch.kevin.simulators.nz;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class NZ_CompModels {
	
	private static final String SUB_SECTS_PATH = "/simulators/nz_comp_subsects/fault_sections.geojson";
	
	public static enum NZ_CompFaultModel implements RupSetFaultModel {
		MODEL_2023;

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
			return "NZNSHM23";
		}

		@Override
		public String getName() {
			return "NZ NSHM23 Fault Model";
		}

		@Override
		public List<? extends FaultSection> getFaultSections() throws IOException {
			throw new UnsupportedOperationException();
		}

		@Override
		public RupSetDeformationModel getDefaultDeformationModel() {
			return null;
		}
		
	}
	
	public static enum NZ_CompDefModel implements RupSetDeformationModel {
		MODEL_2023;
		
		private List<GeoJSONFaultSection> subSects = null;

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
			return "NZNSHM23";
		}

		@Override
		public String getName() {
			return "NZ NSHM23 Deformation Model";
		}

		@Override
		public boolean isApplicableTo(RupSetFaultModel faultModel) {
			return faultModel instanceof NZ_CompFaultModel;
		}

		@Override
		public synchronized List<? extends FaultSection> build(RupSetFaultModel faultModel) throws IOException {
			Preconditions.checkState(faultModel == NZ_CompFaultModel.MODEL_2023);
			if (subSects == null) {
				InputStream is = NZ_CompModels.class.getResourceAsStream(SUB_SECTS_PATH);
				BufferedReader bRead = new BufferedReader(new InputStreamReader(is));
				subSects = GeoJSONFaultReader.readFaultSections(bRead);
			}
			return subSects;
		}
		
	}
	
	public static void main(String[] args) throws IOException {
		List<? extends FaultSection> subSects = NZ_CompDefModel.MODEL_2023.build(NZ_CompFaultModel.MODEL_2023);
		for (FaultSection sect : subSects)
			System.out.println(sect.getSectionId()+". "+sect.getSectionName());
	}

}
