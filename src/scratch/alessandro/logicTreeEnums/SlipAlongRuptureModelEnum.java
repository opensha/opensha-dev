package scratch.alessandro.logicTreeEnums;

public enum SlipAlongRuptureModelEnum {
	
	UNIFORM(	"Uniform",			"Uni"),	// "Uniform/Boxcar (Dsr=Dr)"
	TAPERED(	"Tapered Ends",		"Tap");
	
	private String name, shortName;
	
	private SlipAlongRuptureModelEnum(String name, String shortName) {
		this.name = name;
		this.shortName = shortName;
	}
	
	public String getName() {
		return name;
	}
	
	public String getShortName() {
		return shortName;
	}

	public String getBranchLevelName() {
		return "Slip Along Rupture Model (Dsr)";
	}
}