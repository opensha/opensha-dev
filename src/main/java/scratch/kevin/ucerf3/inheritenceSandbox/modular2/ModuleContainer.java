package scratch.kevin.ucerf3.inheritenceSandbox.modular2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import org.apache.commons.lang3.ClassUtils;
import org.opensha.commons.data.Named;

import com.google.common.base.Preconditions;

/**
 * Container of {@link OpenSHA_Module} instances. Modules are added via {@link #addModule(OpenSHA_Module)}
 * and retrieved by their classes via {@link #getModule(Class)}.
 * <p>
 * When you add a module, all super-classes and super-interfaces of your module that also implement
 * {@link OpenSHA_Module} will be registered, so if you have the following classes:
 * 
 * <pre>
 * abstract class AbstractCustomModule implements OpenSHA_Module {
 * 	// abstract module code
 * }
 * 
 *  class CustomModuleImpl extends AbstractCustomModule {
 * 	// concrete module code
 * }
 * </pre>
 * 
 * ...and you add a module that is an instance of {@code CustomModuleImpl}, you can retrieve it via either
 * {@code getModule(CustomModuleImpl.class)} or {@code getModule(AbstractCustomModule.class)}. Helper classes or
 * interfaces that should never be mapped to a concrete module implementation should be marked with the
 * {@link ModuleHelper} annotation, and will be excluded from any mappings. 
 * <p>
 * TODO: consider synchronization
 * 
 * @author kevin
 *
 * @param <E>
 */
public class ModuleContainer<E extends OpenSHA_Module> {
	
	protected List<E> modules;
	private Map<Class<? extends E>, E> mappings;
	
	private List<Callable<E>> availableModules;
	private Map<Class<? extends E>, Callable<E>> availableMappings;
	
	public ModuleContainer() {
		modules = new ArrayList<>();
		mappings = new HashMap<>();
		availableModules = new ArrayList<>();
		availableMappings = new HashMap<>();
	}
	
	/**
	 * Helper method to determine if a module exists that maps to the given class, equivalent to
	 * {@code getModule(clazz) != null}.
	 * 
	 * @param clazz
	 * @return true if a module exists for the given class, false otherwise
	 */
	public boolean hasModule(Class<? extends E> clazz) {
		return getModule(clazz) != null;
	}
	
	/**
	 * Retrieves a module that matches the given class if it exists, otherwise null.
	 * 
	 * @param <M> the type to be returned
	 * @param clazz module class to look up
	 * @return module mapping the given class, or null if no match
	 */
	@SuppressWarnings("unchecked")
	public <M extends E> M getModule(Class<M> clazz) {
		E module = mappings.get(clazz);
		if (module == null && !availableMappings.isEmpty()) {
			// see if we have it, and then load it lazily
			Callable<E> call = availableMappings.get(clazz);
			if (call != null && loadAvilableModule(call))
				return getModule(clazz);
		}
		return (M)module;
	}
	
	/**
	 * Adds the given module to this container, and maps any applicable super-classes/interfaces (that also implement
	 * {@link OpenSHA_Module} and aren't tagged with the {@link ModuleHelper} annotation) to this instance.
	 * 
	 * @param module
	 */
	public void addModule(E module) {
		// check for duplicates
		for (int m=modules.size(); --m>=0;) {
			E oModule = modules.get(m);
			if (oModule.getClass().equals(module.getClass())) {
				debug("Overriding previous modlue: "+oModule.getName());
				removeMappings(modules.remove(m));
			}
		}
		
		modules.add(module);
		mapModule(module, module.getClass());
		
		for (Class<?> clazz : ClassUtils.getAllSuperclasses(module.getClass()))
			if (isValidModuleSubclass(clazz))
				// this is a super-class
				mapModule(module, clazz);
		for (Class<?> clazz : ClassUtils.getAllInterfaces(module.getClass()))
			if (isValidModuleSubclass(clazz))
				// this is a super-interface
				mapModule(module, clazz);
	}
	
	/**
	 * 
	 * @param clazz
	 * @return true if this is an OpenSHA_Module instance and is not a ModuleHelper
	 */
	private static boolean isValidModuleSubclass(Class<?> clazz) {
		if (!OpenSHA_Module.class.isAssignableFrom(clazz))
			// this is not a module
			return false;
		if (clazz.getAnnotation(ModuleHelper.class) != null)
			// this is a helper class that should not be mapped, skip
			return false;
		return true;
	}
	
	/**
	 * Maps the given module as an implementation of the given class
	 * 
	 * @param module
	 * @param clazz
	 */
	@SuppressWarnings("unchecked")
	private void mapModule(E module, Class<?> clazz) {
		Preconditions.checkState(clazz.getAnnotation(ModuleHelper.class) == null,
				"Cannot map a class that implements @ModuleHelper: %s", clazz.getName());
		if (mappings.containsKey(clazz))
			debug("Overriding module type '"+clazz.getName()+"' with: "+module.getName());
		else
			debug("Mapping module type '"+clazz.getName()+"' to: "+module.getName());
		Preconditions.checkState(OpenSHA_Module.class.isAssignableFrom(clazz));
		mappings.put((Class<E>)clazz, module);
	}
	
	/**
	 * Remove the given module and any mappings to that module.
	 * 
	 * @param module
	 * @return true if the module was present
	 */
	public boolean removeModule(E module) {
		boolean ret = modules.remove(module);
		if (ret)
			removeMappings(module);
		return ret;
	}
	
	/**
	 * Remove any mappings to the given module class (or its subclasses). Also removes any available modules that
	 * map to this class.
	 * 
	 * @param clazz
	 * @return true if any mappings were removed
	 */
	public boolean removeModuleInstances(Class<? extends E> clazz) {
		E module = mappings.get(clazz);
		boolean ret = false;
		if (module != null) 
			ret = removeModule(module);
		
		if (removeAvailableModuleInstance(clazz))
			ret = true;
		
		return ret;
	}
	
	/**
	 * Removes all modules (including any available modules not yet loaded)
	 */
	public void clearModules() {
		modules.clear();
		mappings.clear();
		availableMappings.clear();
		availableModules.clear();
	}
	
	private void removeMappings(E module) {
		List<Class<? extends E>> oldMappings = new ArrayList<>();
		for (Class<? extends E> clazz : mappings.keySet())
			if (mappings.get(clazz).equals(module))
				oldMappings.add(clazz);
		for (Class<? extends E> oldMapping : oldMappings)
			mappings.remove(oldMapping);
	}
	
	/*
	 * Available modules that can be lazily initialized
	 */
	
	/**
	 * Adds an available module that will be lazily loaded when {@link ModuleContainer#getModule(Class)} or
	 * {@link ModuleContainer#hasModule(Class)} is called and no module is presently loaded that matches the given class.
	 * 
	 * @param call
	 * @param moduleClass
	 */
	public void addAvailableModule(Callable<E> call, Class<E> moduleClass) {
		availableModules.add(call);
		mapAvailableModule(call, moduleClass);
		
		for (Class<?> clazz : ClassUtils.getAllSuperclasses(moduleClass))
			if (isValidModuleSubclass(clazz))
				// this is a super-class
				mapAvailableModule(call, clazz);
		for (Class<?> clazz : ClassUtils.getAllInterfaces(moduleClass))
			if (isValidModuleSubclass(clazz))
				// this is a super-interface
				mapAvailableModule(call, clazz);
	}
	
	/**
	 * Loads and maps any available modules that have not yet been initialized
	 * 
	 * @see {@link #addAvailableModule(Callable, Class)}
	 */
	public void loadAllAvailableModules() {
		// wrap in new list, as the load method modifies this list
		List<Callable<E>> available = new ArrayList<>(availableModules);
		for (Callable<E> call : available)
			loadAvilableModule(call);
	}
	
	/**
	 * Loads the given available module
	 * 
	 * @param call (must have already been registered via {@linkp #addAvailableModule(Callable, Class)})
	 * @return true if loading succeeded
	 * @throws IllegalStateException if call is not already registered as an available module
	 */
	public boolean loadAvilableModule(Callable<E> call) {
		Preconditions.checkState(availableModules.remove(call));
		E module = null;
		try {
			debug("Lazily loading available module...");
			module = call.call();
		} catch (Exception e) {
			e.printStackTrace();
			debug("WARNING: failed to lazily load a module (see exception above)", true);
		}
		
		// remove this available module, whether or not it was successful
		availableModules.remove(call);
		List<Class<? extends E>> oldMappings = new ArrayList<>();
		for (Class<? extends E> oClazz : mappings.keySet())
			if (availableMappings.get(oClazz).equals(call))
				oldMappings.add(oClazz);
		for (Class<? extends E> oldMapping : oldMappings)
			availableMappings.remove(oldMapping);
		
		if (module != null) {
			// register it and report success
			addModule(module);
			return true;
		}
		// failed
		return false;
	}
	
	@SuppressWarnings("unchecked")
	private void mapAvailableModule(Callable<E> call, Class<?> clazz) {
		Preconditions.checkState(clazz.getAnnotation(ModuleHelper.class) == null,
				"Cannot map a class that implements @ModuleHelper: %s", clazz.getName());
		if (availableMappings.containsKey(clazz))
			debug("Overriding available module with type: "+clazz.getName());
		else
			debug("Mapping available module with type: "+clazz.getName());
		Preconditions.checkState(OpenSHA_Module.class.isAssignableFrom(clazz));
		availableMappings.put((Class<E>)clazz, call);
	}
	
	/**
	 * Removes any available modules that have been added for the given class
	 * 
	 * @param clazz
	 * @return true if an available module was removed
	 */
	public boolean removeAvailableModuleInstance(Class<? extends E> clazz) {
		Callable<E> call = availableMappings.get(clazz);
		if (call != null) {
			availableModules.remove(call);
			removeAvailableMappings(call);
			return true;
		}
		return false;
	}
	
	private void removeAvailableMappings(Callable<E> call) {
		List<Class<? extends E>> oldMappings = new ArrayList<>();
		for (Class<? extends E> clazz : mappings.keySet())
			if (availableMappings.get(clazz).equals(call))
				oldMappings.add(clazz);
		for (Class<? extends E> oldMapping : oldMappings)
			availableMappings.remove(oldMapping);
	}
	
	/**
	 * @return unmodifiable view of the current modules (not including any available modules not-yet loaded)
	 */
	public List<E> getModules() {
		return getModules(false);
	}
	
	/**
	 * @param loadAvailable if true, any available but not-yet loaded modules will be loaded first
	 * @return unmodifiable view of the current modules
	 */
	public List<E> getModules(boolean loadAvailable) {
		if (loadAvailable)
			loadAllAvailableModules();
		return Collections.unmodifiableList(modules);
	}
	
	/**
	 * 
	 * @return unmodifiable view of the current available module loaders
	 */
	public List<Callable<E>> getAvailableModules() {
		return Collections.unmodifiableList(availableModules);
	}
	
	/**
	 * This returns a unique prefix to be used if this container also a {@link OpenSHA_Module} and is written as member
	 * of a {@link ModuleArchive}. This allows nested file structures within an archive. Default implementation returns
	 * null, and must be overridden to supply a non-empty prefix if this is ever included as a module within a parent
	 * archive.
	 * <p>
	 * Implementations may wish to return a string ending with a forward slash ('/'), which will create a new directory
	 * for all files within an archive.
	 * 
	 * @return prefix that will be applied to files for all modules if this container is written to an archive
	 */
	protected String getNestingPrefix() {
		return null;
	}
	
	private void debug(String message) {
		debug(message, false);
	}
	
	private void debug(String message, boolean err) {
		if (this instanceof Named)
			message = ((Named)this).getName()+":\t"+message;
		if (err)
			System.err.println(message);
		else
			System.out.println(message);
	}

}
