package scratch.ned.nshm23;

/**
 * This class stores a value for N and T, and optionally T2 and/or a weight (which default
 * to NaN and 1.0, respectively).  If T2 is provided, T is assumed to be uncertain with a 
 * uniform probability between T and T2. The weight parameter can be used to store a logic-
 * tree branch weight.  See discussion above on how N can be interpreted (either as number of
 * events or number of intervals).
 * @author field
 *
 */
class NinT_Data {
	final int n;
	final double t;
	final double t2;
	final double weight;
	
	public NinT_Data(int n, double t) {
		this(1.0,n,t,Double.NaN);
	}
	
	public NinT_Data(int n, double t, double t2) {
		this(1.0,n,t,t2);
	}
	
	public NinT_Data(double weight, int n, double t) {
		this(weight,n,t,Double.NaN);
	}
	
	public NinT_Data(double weight, int n, double t, double t2) {
		this.weight=weight;
		this.n=n;
		this.t=t;
		this.t2=t2;
	}
}
