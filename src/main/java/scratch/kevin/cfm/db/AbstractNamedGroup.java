package scratch.kevin.cfm.db;

import org.opensha.commons.data.ShortNamed;

abstract class AbstractNamedGroup implements ShortNamed {
	
	private int id;
	private String name;
	private String shortName;
	
	public AbstractNamedGroup(int id, String name, String shortName) {
		this.id = id;
		this.name = name;
		this.shortName = shortName;
	}

	@Override
	public String getName() {
		return name;
	}
	@Override
	public String getShortName() {
		return shortName;
	}
	public int getID() {
		return id;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + id;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result
				+ ((shortName == null) ? 0 : shortName.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AbstractNamedGroup other = (AbstractNamedGroup) obj;
		if (id != other.id)
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (shortName == null) {
			if (other.shortName != null)
				return false;
		} else if (!shortName.equals(other.shortName))
			return false;
		return true;
	}
	
	// TODO add methods to get from DB
	

}
