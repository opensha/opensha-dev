package scratch.kevin.tdProbModelPlayground;

import java.util.EnumSet;

import org.opensha.commons.data.WeightedList;

/**
 * Enum of probability model options, similar to current ProbabilityModelOptions
 */
public enum FSS_ERF_ProbabilityModels {

	POISSON("Poisson") {
		@Override
		public FSS_ERF_ProbabilityModel.Poisson getProbabilityModel() {
			return new FSS_ERF_ProbabilityModel.Poisson();
		}
	},
	UCERF3_BPT("UCERF3 BPT") {
		@Override
		public UCERF3_ProbabilityModel getProbabilityModel() {
			return new UCERF3_ProbabilityModel(
					// initialize with U3 middle aperiodicity and allow only the 3 UCERF3 aperiodicity branches
					AperiodicityModels.UCERF3_MIDDLE, AperiodicityModels.UCERF3_MODELS,
					// allow only BPT
					RenewalDistributions.BPT, EnumSet.of(RenewalDistributions.BPT));
		}
	},
	NSHM26("NSHM (2026)") {
		@Override
		public UCERF3_ProbabilityModel getProbabilityModel() {
			return new UCERF3_ProbabilityModel(
					// initialize with NSHM26 middle aperiodicity and allow only the NSHM26 aperiodicity branches
					AperiodicityModels.NSHM26_MIDDLE, AperiodicityModels.NSHM26_MODELS,
					// initialize with BPT but allow any of the renewal model distributions
					RenewalDistributions.BPT, EnumSet.allOf(RenewalDistributions.class));
		}
	},
	UCERF3_PREF_BLEND("UCERF3 Preferred Blend") {
		@Override
		public FSS_ERF_ProbabilityModel getProbabilityModel() {
			WeightedList<FSS_ERF_ProbabilityModel> models = new WeightedList<>(4);
			// hardcode these to not allow changing the aperiodicity or renewal models
			models.add(new UCERF3_ProbabilityModel(
					AperiodicityModels.UCERF3_LOW, EnumSet.of(AperiodicityModels.UCERF3_LOW),
					RenewalDistributions.BPT, EnumSet.of(RenewalDistributions.BPT)), 0.1);
			models.add(new UCERF3_ProbabilityModel(
					AperiodicityModels.UCERF3_MIDDLE, EnumSet.of(AperiodicityModels.UCERF3_MIDDLE),
					RenewalDistributions.BPT, EnumSet.of(RenewalDistributions.BPT)), 0.4);
			models.add(new UCERF3_ProbabilityModel(
					AperiodicityModels.UCERF3_HIGH, EnumSet.of(AperiodicityModels.UCERF3_HIGH),
					RenewalDistributions.BPT, EnumSet.of(RenewalDistributions.BPT)), 0.3);
			models.add(new FSS_ERF_ProbabilityModel.Poisson(), 0.2);
			return new FSS_ERF_ProbabilityModel.WeightedCombination(models);
		}
	},
	WG02("WGCEP (2002)") {
		@Override
		public WG02_ProbabilityModel getProbabilityModel() {
			return new WG02_ProbabilityModel();
		}
	};
	
	private String name;

	private FSS_ERF_ProbabilityModels(String name) {
		this.name = name;
	}
	
	@Override
	public String toString() {
		return name;
	}
	
	public abstract FSS_ERF_ProbabilityModel getProbabilityModel();
}
