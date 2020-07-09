package scratch.kevin.ucerf3.eal.branches;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

public enum U3_EAL_GM_Variability implements LogicTreeBranchNode<U3_EAL_GM_Variability> {
	NEG2("-2", "Neg2", 0.1, -2d),
	NEG1("-1", "Neg1", 0.1, -1d),
	ZERO("0", "Zero", 0.6, 0d),
	POS1("+1", "Pos1", 0.1, 1d),
	POS2("+2", "Pos2", 0.1, 2d);

	private String name;
	private String shortName;
	private double weight;
	private double val;
	
	private U3_EAL_GM_Variability(String name, String shortName, double weight, double val) {
		this.name = name;
		this.shortName = shortName;
		this.weight = weight;
		this.val = val;
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
	public double getRelativeWeight(InversionModels im) {
		return weight;
	}

	@Override
	public String encodeChoiceString() {
		return getShortName();
	}

	@Override
	public String getBranchLevelName() {
		return "Combined GM Variability";
	}
	
	public double calcGMVarLoss(double lambda) {
		return calcGMVarLoss(lambda, val);
	}
	
	// constants for GM variability
	private static final double a_sub_L = 0.1992;
	private static final double b_sub_L = 1.0696;
	private static final double c_sub_L = 3.1055;
	private static final double d_sub_L = -0.0935;
	private static final double e_sub_L = 0.091;

	private static double calcGMVarLoss(double lambda, double epsilon_sub_L) {
		if (lambda == 0d)
			return 0d;
		// eqn 35
		double theta_sub_L = a_sub_L * Math.pow(lambda, b_sub_L);
		// eqn 37
		double beta_sub_L0 = c_sub_L * Math.pow(theta_sub_L, d_sub_L);
		// eqn 38
		double beta_sub_L1 = e_sub_L;
		// eqn 36
		double beta_sub_L = Math.sqrt(beta_sub_L0*beta_sub_L0 + beta_sub_L1*beta_sub_L1);
		// eqn 34
		double A_sub_L = (1d/lambda) * theta_sub_L * Math.exp(epsilon_sub_L * beta_sub_L);
		// eqn 33
		return A_sub_L * lambda;
	}
	
	public static void main(String[] args) {
		double lambda = 10000;
		
		System.out.println("Nominal loss: "+lambda);
		for (U3_EAL_GM_Variability var : values()) {
			double varLoss = var.calcGMVarLoss(lambda);
			double ratio = varLoss/lambda;
			System.out.println(var.name+": varLoss = "+(float)varLoss+",\tratio="+(float)ratio);
		}
	}

}
