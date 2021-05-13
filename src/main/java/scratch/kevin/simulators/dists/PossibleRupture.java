package scratch.kevin.simulators.dists;

public class PossibleRupture {
	
	private double prob;
	private double eventTimeYears;
	
	public PossibleRupture(double prob, double eventTimeYears) {
		super();
		this.prob = prob;
		this.eventTimeYears = eventTimeYears;
	}

	public double getProb() {
		return prob;
	}

	public double getEventTimeYears() {
		return eventTimeYears;
	}

}
