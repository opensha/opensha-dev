package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;


/**
 * Operation payload for writing message to console.
 * Author: Michael Barall 05/05/2018.
 */
public class OpConsoleMessage extends DBPayload {

	//----- Constants and variables -----

	// The message to be written to the console.
	// Cannot be null, but can be an empty string.

	public String message;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpConsoleMessage () {}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 26001;

	private static final String M_VERSION_NAME = "OpConsoleMessage";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalString ("message", message);

		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Superclass

		super.do_umarshal (reader);

		// Contents

		message = reader.unmarshalString ("message");

		return;
	}

	// Marshal object.

	@Override
	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	@Override
	public OpConsoleMessage unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

}
