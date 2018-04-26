package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalWriter;

import org.bson.types.ObjectId;

/**
 * Holder for database record payload.
 * Author: Michael Barall 04/08/2018.
 *
 * In the AAFS server, some database records hold a payload, which the user can
 * access by marshaling and unmarshaling.  This object holds a payload in the
 * form that it is stored in the database.  It can be used to transfer a payload
 * from one record to another.  It implements MarshalWriter so it can be passed
 * in to the same routines that are used to supply marshaled data.

 *
 * Only code very close to the database engine should create these objects or
 * access their contents.  All other code should treat these objects as opaque.
 * (If this was C++, everything would be private and the classes that need access
 * would be friends, but Java does not have an analogous mechanism.)
 */
public class RecordPayload implements MarshalWriter {

	// The JSON string.

    private String json_string;

	// Constructor saves the JSON string.

	public RecordPayload (String json_string) {
		this.json_string = json_string;
	}

	// Get the JSON string.

	public String get_json_string () {
		return json_string;
	}




	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "RecordPayload: " + ((json_string == null) ? ("null") : ("len = " + json_string.length()));
		return str;
	}




	//----- Empty Implementation of MarshalWriter


	/**
	 * Begin a map context.
	 */
	@Override
	public void marshalMapBegin (String name) {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}

	/**
	 * End a map context.
	 */
	@Override
	public void marshalMapEnd () {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}

	/**
	 * Begin an array context, specify the array size.
	 */
	@Override
	public void marshalArrayBegin (String name, int array_size) {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}

	/**
	 * End an array context.
	 */
	@Override
	public void marshalArrayEnd () {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}

	/**
	 * Marshal a long.
	 */
	@Override
	public void marshalLong (String name, long x) {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}

	/**
	 * Marshal a double.
	 */
	@Override
	public void marshalDouble (String name, double x) {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}

	/**
	 * Marshal a string.  (Null strings are not allowed.)
	 */
	@Override
	public void marshalString (String name, String x) {
		throw new UnsupportedOperationException ("RecordPayload does not support writing");
	}



}
