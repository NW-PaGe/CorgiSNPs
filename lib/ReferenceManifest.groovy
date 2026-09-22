//
// Reference manifest loader.
//
// Place this file in the pipeline's lib/ directory; Nextflow loads it automatically.
//
//   def manifest = ReferenceManifest.load(params.reference_db)          // validate
//   def manifest = ReferenceManifest.load(params.reference_db, false)   // skip validation
//
//   manifest.records  one map per subtype (see below)
//   manifest.errors   list of problems; empty when the directory is valid
//
// The directory holds manifest.yml plus one sub-directory per species entry,
// named by the entry's 'name', with one folder per file field:
//
//   reference_db/
//     manifest.yml
//     candidozyma_auris/
//       assembly/GCA_016772135.1.fna.gz
//       annotation/GCA_016772135.1.gff.gz
//
// File fields (assembly, annotation) are looked up first as written (absolute, or
// relative to the launch directory), then at <reference_db>/<name>/<field>/<file>.
//
// Each record contains every field from its species entry (except 'subtypes')
// plus every field from the subtype itself; if both define a field, the
// subtype's value wins. 'species' and 'subtype' are always lists, file fields
// are resolved paths, and 'reference' is the resolved assembly.
//
// Required: 'name', 'species' and 'subtypes' on each entry; 'subtype' and
// 'assembly' on each subtype. Everything else is optional and passed through.
//
// AMR: a subtype with an 'amr' block must also have an 'annotation'. Each
// record gets a boolean 'primary': the subtype marked 'primary: true', or, if
// none is marked, the only subtype of that species with 'amr'. Marking more
// than one, marking one without 'amr', or leaving several 'amr' subtypes with
// none marked is an error.
//
// When validation fails, records is empty and errors lists every problem found.
//

import java.nio.file.Path
import groovy.yaml.YamlSlurper
import nextflow.Nextflow

class ReferenceManifest {

    static final String MANIFEST_NAME = 'manifest.yml'

    // Subtype fields that name files; these are resolved to paths
    static final List<String> FILE_FIELDS = [ 'assembly', 'annotation' ]

    // Allowed characters for an entry's name, since it is used as a directory name
    static final String NAME_PATTERN = '^[A-Za-z0-9._-]+$'

    static Map load(dir, boolean validate = true) {
        Path root = Nextflow.file(dir.toString()) as Path
        Path manifest = root.resolve(MANIFEST_NAME)
        if( !manifest.exists() )
            return result([], [ "No ${MANIFEST_NAME} found in ${dir}" ])

        def entries
        try {
            entries = new YamlSlurper().parse(manifest)
        }
        catch( Exception e ) {
            return result([], [ "Could not parse ${manifest}: ${e.message}" ])
        }
        if( !(entries instanceof List) )
            return result([], [ "${manifest} must contain a list of entries" ])

        List<Map> records = []
        List errors = []
        Map seenNames = [:]

        entries.eachWithIndex { entry, i ->
            if( !(entry instanceof Map) ) {
                errors << "Entry ${i + 1}: expected a mapping with 'name', 'species' and 'subtypes' keys"
                return
            }

            def name = entry.name?.toString()?.trim()
            def species = asStringList(entry.species)
            def label = "Entry ${i + 1}" + (name ? " (${name})" : (species ? " (${species[0]})" : ''))
            def subtypes = entry.subtypes instanceof List ? entry.subtypes : []

            if( validate ) {
                if( !name )
                    errors << "${label}: 'name' must be set"
                else if( !(name ==~ NAME_PATTERN) )
                    errors << "${label}: 'name' can only contain letters, numbers, '.', '_' and '-'"
                else if( seenNames.containsKey(name) )
                    errors << "${label}: duplicate name '${name}' (first used in entry ${seenNames[name]})"
                else
                    seenNames[name] = i + 1

                if( !species )
                    errors << "${label}: 'species' must be set"
                if( !subtypes )
                    errors << "${label}: 'subtypes' must be a non-empty list"
            }

            Path speciesDir = name ? root.resolve(name) : root

            // Subtype name -> position where it was first used, for the uniqueness check
            Map seen = [:]
            int first = records.size()

            subtypes.eachWithIndex { st, j ->
                def stLabel = "${label}, subtype ${j + 1}"
                if( !(st instanceof Map) ) {
                    errors << "${stLabel}: expected a mapping with 'subtype' and 'assembly' keys"
                    return
                }

                def names = asStringList(st.subtype)
                if( validate ) {
                    if( !names )
                        errors << "${stLabel}: 'subtype' must be set"
                    names.each { n ->
                        if( seen.containsKey(n) )
                            errors << "${stLabel}: duplicate subtype '${n}' (first used in subtype ${seen[n]})"
                        else
                            seen[n] = j + 1
                    }
                }

                // Resolve file fields to paths
                Map files = [:]
                FILE_FIELDS.each { field ->
                    def value = st[field]?.toString()?.trim()
                    if( !value ) {
                        if( validate && field == 'assembly' )
                            errors << "${stLabel}: 'assembly' must be set"
                        return
                    }
                    def found = candidates(value, speciesDir, field).find { it.exists() }
                    if( validate && !found )
                        errors << "${stLabel}: ${field} not found. Tried:\n" +
                            candidates(value, speciesDir, field).collect { "      ${it}" }.join('\n')
                    files[field] = found
                }

                def shared = entry.findAll { k, v -> k != 'subtypes' }
                records << shared + st + files + [ name: name, species: species, subtype: names, reference: files.assembly ]
            }

            // AMR checks and primary selection for this species
            List entryRecords = records[first..<records.size()]
            List withAmr = entryRecords.findAll { r -> hasAmr(r) }
            List flagged = entryRecords.findAll { r -> isTrue(r.primary) }

            if( validate ) {
                withAmr.findAll { r -> !r.annotation }.each { r ->
                    errors << "${label}, subtype '${r.subtype[0]}': 'amr' requires an 'annotation'"
                }
                withAmr.each { r ->
                    def targets = r.amr instanceof List ? r.amr : [ r.amr ]
                    targets.eachWithIndex { t, n ->
                        if( !(t instanceof Map) || !t.gene?.toString()?.trim() )
                            errors << "${label}, subtype '${r.subtype[0]}': amr target ${n + 1} must set 'gene'"
                    }
                }
                if( flagged.size() > 1 )
                    errors << "${label}: only one subtype can be primary, found ${flagged.collect { r -> r.subtype[0] }.join(', ')}"
                else if( flagged && !hasAmr(flagged[0]) )
                    errors << "${label}: primary subtype '${flagged[0].subtype[0]}' must have 'amr' and 'annotation'"
                else if( !flagged && withAmr.size() > 1 )
                    errors << "${label}: ${withAmr.size()} subtypes have 'amr' (${withAmr.collect { r -> r.subtype[0] }.join(', ')}); mark one with 'primary: true'"
            }

            def primary = flagged.size() == 1 ? flagged[0] : (!flagged && withAmr.size() == 1 ? withAmr[0] : null)
            for( int k = first; k < records.size(); k++ )
                records[k] = records[k] + [ primary: records[k].is(primary) ]
        }

        return errors ? result([], errors) : result(records, [])
    }

    //
    // Locations checked for a file field, in order:
    //   1. the path exactly as written (absolute, or relative to the launch directory)
    //   2. <reference_db>/<name>/<field>/<path>, e.g. candidozyma_auris/annotation/GCA_016772135.1.gff.gz
    //
    static List<Path> candidates(String value, Path speciesDir, String field) {
        return [ Nextflow.file(value) as Path, speciesDir.resolve(field).resolve(value) ]
    }

    static boolean hasAmr(Map record) {
        def amr = record.amr
        return amr instanceof Collection ? !amr.isEmpty() : amr != null
    }

    static List<String> amrGenes(Map record) {
        if( !hasAmr(record) )
            return []
        def targets = record.amr instanceof List ? record.amr : [ record.amr ]
        return targets
            .findAll { t -> t instanceof Map }
            .collect { t -> t.gene?.toString()?.trim() }
            .findAll { g -> g }
            .unique()
    }

    static boolean isTrue(value) {
        return value instanceof Boolean ? value : value?.toString()?.trim()?.toLowerCase() in [ 'true', 'yes' ]
    }

    static List<String> asStringList(value) {
        def items = value instanceof List ? value : (value == null ? [] : [value])
        return items.findAll { it != null }.collect { it.toString().trim() }.findAll { it }
    }

    static Map result(List records, List errors) {
        return [ records: records, errors: errors.collect { it.toString() } ]
    }
}