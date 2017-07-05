package scratch.kevin.markov;

import java.util.Arrays;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;

public class IndicesKey {
	private static Joiner j = Joiner.on(",");
	private int hashCode;
	int[] indices;
	public IndicesKey(int[] indices) {
		this.indices = indices;
		this.hashCode = Arrays.hashCode(indices);
	}
	@Override
	public int hashCode() {
		return hashCode;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		IndicesKey other = (IndicesKey) obj;
		return Arrays.equals(indices, other.indices);
	}
	
	public String toString() {
		return "IndicesKey["+j.join(Ints.asList(indices))+"]";
	}
	
	public int[] getIndices() {
		return indices;
	}
}