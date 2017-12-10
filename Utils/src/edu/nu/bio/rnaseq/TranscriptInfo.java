package edu.nu.bio.rnaseq;

class TranscriptInfo {
	String transcriptId;
	String geneID;
	Integer cluster;

	public void update(TranscriptInfo info) {
		if (this.transcriptId == null) {
			this.transcriptId = info.transcriptId;
		}
		if (this.geneID == null) {
			this.geneID = info.geneID;
		}
		if (this.cluster == null) {
			this.cluster = info.cluster;
		}
	}

	@Override
	public boolean equals(Object obj) {

		if (((TranscriptInfo) obj).cluster == this.cluster
				&& ((TranscriptInfo) obj).transcriptId
						.equals(this.transcriptId)
				&& ((TranscriptInfo) obj).geneID.equals(this.geneID)) {
			return true;
		} else {
			return false;
		}

	}
	
	
	
	
}