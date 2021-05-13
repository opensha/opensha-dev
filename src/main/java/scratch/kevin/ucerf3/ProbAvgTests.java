package scratch.kevin.ucerf3;

public class ProbAvgTests {

	public static void main(String[] args) {
		int num = 1440*4;
		double duration = 5;
		
		double sumProb = 0;
		double sumRate = 0;
		for (int i=0; i<num; i++) {
			double rate = 1e-2*Math.random();
			double prob = 1d-Math.exp(-rate*duration);
			sumRate += rate;
			sumProb += prob;
		}
		
		double avgProb = sumProb/(double)num;
		double avgRate = sumRate/(double)num;
		double calcAvgProb = 1d-Math.exp(-avgRate*duration);
		
		System.out.println("Avg Prob: "+avgProb);
		System.out.println("Avg Rate: "+avgRate);
		System.out.println("Calc Avg Prob: "+calcAvgProb);
	}

}
