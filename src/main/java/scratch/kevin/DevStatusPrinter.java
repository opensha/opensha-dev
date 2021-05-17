package scratch.kevin;

import java.util.Set;

import org.opensha.commons.util.DevStatus;
import org.opensha.sha.earthquake.ERF_Ref;
import org.opensha.sha.imr.AttenRelRef;

public class DevStatusPrinter {
	
	private static void printERFs() {
		System.out.println("'''ERFs'''\n");
		printStatus(ERF_Ref.get(true, DevStatus.PRODUCTION), DevStatus.PRODUCTION);
		printStatus(ERF_Ref.get(true, DevStatus.DEVELOPMENT), DevStatus.DEVELOPMENT);
		printStatus(ERF_Ref.get(true, DevStatus.EXPERIMENTAL), DevStatus.EXPERIMENTAL);
		printStatus(ERF_Ref.get(true, DevStatus.DEPRECATED), DevStatus.DEPRECATED);
	}
	
	private static void printIMRs() {
		System.out.println("'''IMRs'''\n");
		printStatus(AttenRelRef.get(DevStatus.PRODUCTION), DevStatus.PRODUCTION);
		printStatus(AttenRelRef.get(DevStatus.DEVELOPMENT), DevStatus.DEVELOPMENT);
		printStatus(AttenRelRef.get(DevStatus.EXPERIMENTAL), DevStatus.EXPERIMENTAL);
		printStatus(AttenRelRef.get(DevStatus.DEPRECATED), DevStatus.DEPRECATED);
	}
	
	private static void printStatus(Set<?> refs, DevStatus ds) {
		switch (ds) {
			case PRODUCTION:
				System.out.println("Production (these show up in all applications)");
				break;
			case DEVELOPMENT:
				System.out.println("Development (these show up in nightly builds or when run from trunk)");
				break;
			case EXPERIMENTAL:
				System.out.println("Experimental (these show up in nightly builds or when run from trunk)");
				break;
			case DEPRECATED:
				System.out.println("Deprecated (these never show up in any apps)");
				break;
		}
		for (Object ref : refs)
			System.out.println("* "+ref.toString());
		System.out.println();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		printERFs();
		printIMRs();
	}

}
