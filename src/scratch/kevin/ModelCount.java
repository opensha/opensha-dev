package scratch.kevin;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.util.DevStatus;
import org.opensha.sha.earthquake.ERF_Ref;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Joiner;

public class ModelCount {

	public static void main(String[] args) {
		List<DevStatus[]> levels = new ArrayList<>();
		levels.add(new DevStatus[] {DevStatus.PRODUCTION});
		levels.add(new DevStatus[] {DevStatus.PRODUCTION, DevStatus.DEVELOPMENT});
		levels.add(new DevStatus[] {DevStatus.PRODUCTION, DevStatus.DEVELOPMENT, DevStatus.EXPERIMENTAL});
		
		for (DevStatus[] status : levels) {
			System.out.println(Joiner.on(", ").join(status));
			System.out.println("\tERFs: "+ERF_Ref.get(true, status).size());
			System.out.println("\tGMPEs: "+AttenRelRef.get(status).size());
		}
	}

}
