package scratch.kevin.matBench;

public abstract class AbstractBenchableMatrix2D implements BenchableMatrix2D {
	
	long multTime = 0;

	@Override
	public final BenchableMatrix2D mult(BenchableMatrix2D mat) {
		long start = System.currentTimeMillis();
		BenchableMatrix2D res = doMult(mat);
		multTime += System.currentTimeMillis() - start;
		return res;
	}
	
	public long getMultTime() {
		return multTime;
	}
	
	public void resetMultTime() {
		multTime = 0;
	}
	
	protected abstract BenchableMatrix2D doMult(BenchableMatrix2D mat);

}
