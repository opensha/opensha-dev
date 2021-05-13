package scratch.kevin.cfm.db;

import org.opensha.commons.data.Named;
import org.opensha.commons.util.ClassUtils;

public class CFMFault implements Named {
	
	private String objectName;
	
	private CFMFaultSystem system;
	private CFMFaultZone zone;
	private CFMFaultSection section;
	private CFMFaultName faultName;
	private CFMFaultStrand strand;
	
	public CFMFault(String objectName, CFMFaultSystem system,
			CFMFaultZone zone, CFMFaultSection section, CFMFaultName faultName,
			CFMFaultStrand strand) {
		this.objectName = objectName;
		this.system = system;
		this.zone = zone;
		this.section = section;
		this.faultName = faultName;
		this.strand = strand;
	}

	@Override
	public String getName() {
		return objectName;
	}
	
	public String toString() {
		return objectName+":\t"+getGroupName(CFMFaultSystem.class, system)
				+"\t"+getGroupName(CFMFaultZone.class, zone)
				+"\t"+getGroupName(CFMFaultSection.class, section)
				+"\t"+getGroupName(CFMFaultName.class, faultName)
				+"\t"+getGroupName(CFMFaultStrand.class, strand);
	}
	
	private <E extends AbstractNamedGroup> String getGroupName(Class<E> clazz, E group) {
		String name;
		if (group == null)
			name = "(null)";
		else
			name = group.getShortName();
		return ClassUtils.getClassNameWithoutPackage(clazz)+":"+name;
	}

}
