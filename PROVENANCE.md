# Provenance and limits

The source study used internal identifiers W006/W007 and the project name Waddle (now Avify Dose). Historical names and historical certification terminology remain in frozen records/code; the revised manuscript and README govern scientific interpretation.

Recorded result JSON, reference transcription, both transport harnesses, the reference extractor and W007 protocol are copied byte-for-byte. They contain tally means/standard errors and run metadata, not complete per-history records or OpenMC statepoint HDF5 files. W006/W007 validation grades remain INCOMPLETE-VALIDATION. A successful checker exit means consistency, not a successful clinical or complete benchmark validation.

The W006 public protocol omits its internal patent-strategy discussion between the paragraph beginning “Why it is nearest art on the mechanism” and section 9. Its scientific setup and gates remain present. Its original hash is recorded in provenance/original-sha256.json; the original unredacted file is not distributed. Public check_w006.py verifies the public protocol hash and separately checks the original protocol identifier stored in results. Public check_w007.py accepts the adapted check_w006.py hash; its numerical calculations and all other original dependency checks are unchanged. This checks this public snapshot, not possession of the complete frozen original protocol.

The manuscript audit uses the active Python interpreter rather than the author's private virtualenv path. ledger.md contains labelled scientific excerpts only. Author-specific fallback scratch paths remain in the unchanged transport harness; set W006_SCRATCH when running it. Historical references to omitted internal sessions are provenance context, not required dependencies for the recorded-result audit.

Original hashes are in provenance/original-sha256.json. Public release hashes are in MANIFEST.sha256. A clean Python standard-library verification succeeds without the external Goorley archive; that optional archive-level comparison is explicitly skipped unless W006_ARCHIVE is supplied. No full transport campaign is claimed to have been rerun for this packaging release.
