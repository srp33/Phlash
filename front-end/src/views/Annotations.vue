<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <h1>Annotations</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p>
          Gene calls that have been made from <i>GeneMarkS</i>, <i>Glimmer3</i>,
          <i>ARAGORN</i>, and <i>PHANOTATE</i> are shown below. To have the most
          accurate annotation of this genome manual annotation is necessary. The
          status for each gene is shown below which are mearly guides and
          suggestions. The status does not actually mean that a gene call is
          correct or incorrect. The parameters used for creating the status for
          each gene can be updated in the 'Settings' tab above.
        </p>
        <hr />
        <p><strong>Status Definitions</strong></p>
        <p>
          <strong style="color: #1b9e77">Green:</strong> Pass, the open reading
          frame covers the coding potential.<br />
          <strong style="color: #d95f02">Orange:</strong> Fail, the open reading
          frame does not cover the coding potential.<br />
          <strong style="color: #7570b3">Purple:</strong> tRNA, this gene
          represents a tRNA.<br />
          <strong>S: </strong>Short, the open reading frame is less than
          {{ short }} base pairs.<br />
          <strong>LLG: </strong>Long leading gap, there is more than a
          {{ gap }} base pair gap between this and the previous gene.<br />
          <strong>LLO: </strong>Long leading overlap, the gene overlaps the
          previous gene by more than {{ overlap }} base pairs.<br />
          <strong>LTG: </strong>Long tailing gap, there is more than a
          {{ gap }} base pair gap between this and the next gene.<br />
          <strong>LTO: </strong>Long tailing overlap, the gene overlaps the next
          gene by more than {{ overlap }} base pairs.<br />
          <strong>STG: </strong>Short tailing gap, the gene is on a different
          strand than the previous gene and has a gap less than
          {{ oppositeGap }} base pairs in length.<br />
          <strong>SLG: </strong>Short leading gap, the gene is on a different
          strand than the next gene and has a gap less than
          {{ oppositeGap }} base pairs in length.<br />
          <strong>R: </strong>The gene has a reasonable length and overlap.
        </p>
        <hr />
        <p><strong>Actions</strong></p>
        <p>
          <strong>Annotate: </strong>When clicked, a function and alternate open
          reading frame can be added.<br />
          <strong>Delete: </strong>When clicked, this gene will be temporarily
          removed.<br />
          <strong>Reinstate: </strong>When clicked, the deleted gene will be
          added.<br />
          <strong>White background: </strong>This CDS has already been updated
          but can still be edited.<br />
          In rare occasions when a needed coding sequence is not shown, the 'Add
          CDS' button may be clicked to add a new custom CDS.<br />
          <button class="btn btn-dark btn-action" @click="showAddCDS = true">
            <strong>&#43; Add CDS</strong>
          </button>
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'GeneMap',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav" id="next-top">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
      <div v-if="blastLoading" class="alert alert-primary alert-dismissible">
        <a href="#" class="close" data-dismiss="alert" aria-label="close"
          >&times;</a
        >
        The Blast results are being interpretted which can take several
        minutes.<br />
        If 'Annotate' is clicked for a CDS that does not currently have any
        data, you will be brought back to this page.<br />
        A notification will appear when finished.
      </div>
      <div
        class="alert alert-secondary"
        style="text-align: center"
        v-if="completedGenes !== dnamaster.length"
      >
        You have
        <strong
          >{{ dnamaster.length - completedGenes }}/{{
            dnamaster.length
          }}</strong
        >
        genes remaining.
      </div>
      <div
        class="alert alert-secondary"
        v-if="completedGenes === dnamaster.length"
      >
        Congratulations! You have annotated every CDS. Click 'Next' to see a map
        of the genome.
      </div>
      <div id="annotations" align="center">
        <div
          class="table table-responsive table-secondary"
          style="overflow-y: auto; max-height: 50em"
        >
          <table class="table table-hover" align="center">
            <thead>
              <tr>
                <th scope="col">ID</th>
                <th scope="col">Left</th>
                <th scope="col">Right</th>
                <th scope="col">Strand</th>
                <th width="12.5em" scope="col">Product</th>
                <th scope="col">Status</th>
                <th scope="col">Action</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(curr, index) in dnamaster" :key="index">
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                >
                  {{ curr.id }}
                </td>
                <td v-else>{{ curr.id }}</td>
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else>{{ curr.start }}</td>
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else>{{ curr.stop }}</td>
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else-if="curr.strand === '-'">Complementary</td>
                <td v-else>Direct</td>
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td
                  v-else-if="curr.function === 'None selected'"
                  style="font-size: 1.2em"
                >
                  {{ curr.function }}
                </td>
                <td
                  v-else-if="curr.function.length < 17"
                  style="font-size: 1.2em"
                >
                  {{ curr.function.substring(1, curr.function.length) }}
                </td>
                <td v-else style="font-size: 1.2em">
                  {{ curr.function.substring(1, 14) }}...
                </td>
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else-if="curr.status === 'tRNA'" style="color: #7570b3">
                  {{ getStatus(index) }}
                </td>
                <td v-else-if="curr.status === 'Pass'" style="color: #1b9e77">
                  {{ getStatus(index) }}
                </td>
                <td v-else-if="curr.status === 'Fail'" style="color: #d95f02">
                  {{ getStatus(index) }}
                </td>
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                >
                  <button
                    class="btn btn-outline-dark btn-sm"
                    style="width: 6.25em"
                    @click="reinstate(index)"
                  >
                    <strong>Reinstate</strong>
                  </button>
                </td>
                <td v-else-if="curr.status === 'tRNA'">
                  <button
                    class="btn btn-outline-dark btn-sm"
                    style="width: 6.25em"
                    @click="deleteTRNA(index)"
                  >
                    <strong>Delete</strong>
                  </button>
                </td>
                <td v-else>
                  <router-link
                    :to="{
                      name: 'CDS',
                      params: {
                        phageID: $route.params.phageID,
                        cdsID: curr.id,
                      },
                    }"
                  >
                    <button
                      class="btn btn-outline-dark btn-sm"
                      style="width: 6.25em"
                      v-if="curr.function[0] === '@'"
                    >
                      <strong>Annotate</strong>
                    </button>
                  </router-link>
                  <router-link
                    :to="{
                      name: 'CDS',
                      params: {
                        phageID: $route.params.phageID,
                        cdsID: curr.id,
                      },
                    }"
                  >
                    <button
                      class="btn btn-dark btn-sm"
                      style="width: 6.25em"
                      v-if="curr.function === 'None selected'"
                    >
                      <strong>Annotate</strong>
                    </button>
                  </router-link>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
      <div class="alert alert-secondary">
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'GeneMap',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav" id="next-top">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
    <b-modal
      v-model="showAddCDS"
      ref="addCDSModal"
      id="addCDS-modal"
      title="Add CDS"
      hide-footer
    >
      <b-form @submit="onSubmitAdd" align="left">
        <b-form-group label="Left:" label-size="lg" label-for="add-start-input">
          <b-form-input
            id="add-start-input"
            type="number"
            v-model="addCDS.start"
            required
            placeholder="Enter left position"
          ></b-form-input>
        </b-form-group>
        <b-form-group label="Right:" label-size="lg" label-for="add-stop-input">
          <b-form-input
            id="add-stop-input"
            type="number"
            v-model="addCDS.stop"
            required
            placeholder="Enter right position"
          ></b-form-input>
        </b-form-group>
        <b-form-group label="Strand:" label-size="lg">
          <b-form-select
            v-model="addCDS.strand"
            required
            :options="strandOptions"
          ></b-form-select>
        </b-form-group>
        <b-form-group
          label="Force Add:"
          label-size="lg"
          label-for="add-force-input"
        >
          <b-form-checkbox
            id="add-force-input"
            type="checkbox"
            v-model="addCDS.force"
          >
            Do not check this unless you are sure you want to add the CDS. This
            could remove another CDS if a duplicate ID is given or create a CDS
            that does not represent an actual gene.
          </b-form-checkbox>
        </b-form-group>
        <hr />
        <b-button type="submit" class="mt-3" block style="margin-top: 0em">
          <strong>Submit</strong>
        </b-button>
      </b-form>
    </b-modal>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Annotations",
  components: {
    Loading,
    Navbar,
  },

  data() {
    return {
      addCDS: {
        id: "",
        start: "",
        stop: "",
        strand: null,
        force: false,
        read: [],
      },
      genbankAnnotations: {
        phageName: "",
        source: "",
        organism: "",
        isolationSource: "",
        labHost: "",
        identifiedBy: "",
        authors: "",
        title: "",
        journal: "",
        country: "USA",
        molType: "genomic DNA",
        notes: "complete genome",
        includeNotes: false,
      },
      strandOptions: [
        { value: null, text: "Please select a direction" },
        { value: "+", text: "+ (Direct)" },
        { value: "-", text: "- (Complementary)" },
      ],
      dnamaster: [],
      currCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
      },
      showAddCDS: false,
      pageLoading: true,
      blastLoading: true,
      completedGenes: 0,
      gap: 10,
      overlap: 10,
      oppositeGap: 50,
      short: 200,
    };
  },

  created() {
    this.getData();
  },

  computed: {
    navUpload: function () {
      return true;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      return true;
    },

    navGeneMap: function () {
      return true;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },
  },

  methods: {
    /**
     * Populates dnamaster with data.
     */
    getData() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/none`
        )
        .then((response) => {
          this.dnamaster = response.data.dnamaster;
          this.gap = response.data.gap;
          this.overlap = response.data.overlap;
          this.oppositeGap = response.data.opposite_gap;
          this.short = response.data.short;
          this.pageLoading = false;
          this.genbankAnnotations.phageName = this.$route.params.phageID;
          for (var i = 0; i < this.dnamaster.length; i += 1) {
            if (this.dnamaster[i].function !== 'None selected')
              this.completedGenes += 1;
          }
          this.parseBlast();
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Calls method on back-end that parses and stores all of the BLAST results.
     */
    parseBlast() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/blast`
        )
        .then((response) => {
          this.blastLoading = false;
          if (response.data === 'success') {
            this.$bvToast.toast(
              `All of the BLAST results have finished being interpretted.`,
              {
                title: 'Finished',
                appendToast: false,
              }
            );
          } else if (response.data === 'error') {
            this.$bvToast.toast(
              `An unknown error occurred. Try removing and reuploading the BLAST files. 
              If you ignore this error not all of your BLAST results will be shown. 
              If this error continues, please contact us by visiting the 'contact' tab.`,
              {
                title: 'Error',
                appendToast: false,
              }
            );
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Changes the function of a deleted gene so that it is now visible.
     */
    reinstate(index) {
      if (this.dnamaster[index].function === 'DELETED') {
        this.completedGenes -= 1;
        this.dnamaster[index].function = 'None selected';
      } else {
        this.dnamaster[index].status = 'tRNA';
      }
      const payload = {
        id: this.dnamaster[index].id,
        start: this.dnamaster[index].start,
        stop: this.dnamaster[index].stop,
        strand: this.dnamaster[index].strand,
        function: this.dnamaster[index].function,
        status: this.dnamaster[index].status,
        frame: this.dnamaster[index].frame,
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${this.dnamaster[index].id}`,
          payload
        )
        .then(() => {
          console.log(this.dnamaster[index].function);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Temporarily deletes a tRNA gene.
     */
    deleteTRNA(index) {
      this.dnamaster[index].status = 'trnaDELETED';
      const payload = {
        id: this.dnamaster[index].id,
        start: this.dnamaster[index].start,
        stop: this.dnamaster[index].stop,
        strand: this.dnamaster[index].strand,
        function: this.dnamaster[index].function,
        status: 'trnaDELETED',
        frame: this.dnamaster[index].frame,
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${this.dnamaster[index].id}`,
          payload
        )
        .then(() => {
          console.log(this.dnamaster[index].function);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Creates the status for a given cds.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getStatus(index) {
      if (this.dnamaster[index].status === 'tRNA') {
        return 'tRNA';
      }
      var status = '';
      if (
        this.dnamaster[index].stop - this.dnamaster[index].start <
        this.short
      ) {
        status += 'S | ';
      }
      if (index === 0) {
        status += this.getTailingOverlap(index);
      } else if (index === this.dnamaster.length - 1) {
        status += this.getLeadingOverlap(index);
      } else {
        status += this.getLeadingOverlap(index);
        status += this.getTailingOverlap(index);
      }
      if (status === '') {
        return 'R';
      }
      return status.substring(0, status.length - 3);
    },

    /**
     * Determines if there is a large gap or overlap between the current gene and the next.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getTailingOverlap(index) {
      var status = '';
      var nextGene = 1;
      var overlap = -10000;
      while (overlap === -10000) {
        if (index + nextGene >= this.dnamaster.length) {
          return status;
        }
        if (this.dnamaster[index + nextGene].function !== 'DELETED') {
          overlap =
            this.dnamaster[index].stop - this.dnamaster[index + nextGene].start;
        } else {
          nextGene += 1;
        }
      }
      if (
        this.dnamaster[index].strand === this.dnamaster[index + nextGene].strand
      ) {
        if (overlap < this.gap * -1) {
          status += 'LTG | ';
        } else if (overlap > this.overlap) {
          status += 'LTO | ';
        }
      } else if (overlap > this.oppositeGap * -1) {
        status += 'STG | ';
      }
      return status;
    },

    /**
     * Determines if there is a large gap or overlap between the current gene and the previous.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getLeadingOverlap(index) {
      var status = '';
      var nextGene = 1;
      var leadingOverlap = -10000;
      while (leadingOverlap === -10000) {
        if (index - nextGene < 0) {
          return status;
        }
        if (this.dnamaster[index - nextGene].function !== 'DELETED') {
          leadingOverlap =
            this.dnamaster[index - nextGene].stop - this.dnamaster[index].start;
        } else {
          nextGene += 1;
        }
      }
      if (
        this.dnamaster[index].strand === this.dnamaster[index - nextGene].strand
      ) {
        if (leadingOverlap < this.gap * -1) {
          status += 'LLG | ';
        } else if (leadingOverlap > this.overlap) {
          status += 'LLO | ';
        }
      } else if (leadingOverlap > this.oppositeGap * -1) {
        status += 'SLG | ';
      }
      return status;
    },

    /**
     * @param {event} evt the submit event on the add cds modal.
     * Adds a new CDS or notifies the user of a failed add.
     */
    onSubmitAdd(evt) {
      evt.preventDefault();
      this.$refs.addCDSModal.hide();
      let read = true;
      const payload = {
        id: 'added',
        start: this.addCDS.start,
        stop: this.addCDS.stop,
        strand: this.addCDS.strand,
        force: this.addCDS.force,
        read, // property shorthand
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/none`,
          payload
        )
        .then((response) => {
          console.log(response.data.message);
          if (response.data.message === 'ID already exists.') {
            this.$bvToast.toast(
              `The CDS already exists. Try again with a different ORF. To ignore this 
            warning and add the CDS, check the 'Force Add' box.`,
              {
                title: 'ADD FAILED',
                autoHideDelay: 5000,
                appendToast: false,
              }
            );
          } else if (response.data.message === 'Not orf.') {
            this.$bvToast.toast(
              `The inputted start and stop locations do not represent an ORF. To ignore this 
            warning and add the CDS, check the 'Force Add' box.`,
              {
                title: 'ADD FAILED',
                autoHideDelay: 5000,
                appendToast: false,
              }
            );
          } else {
            window.location.reload();
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },
  },
};
</script>

<style scoped>
.wrapper {
  margin: 0;
}

h1 {
  margin-top: 0.7em;
}

.nav-btns-wrapper {
  text-align: center;
}

td,
th {
  word-wrap: break-word;
  width: 9.5em;
  font-size: 1.4em;
}

.table-responsive thead th {
  position: sticky;
  top: 0;
  background: white;
  color: black;
  border: darkgray;
  font-size: 1.4em;
}

.btn-nav {
  margin: 0.25em;
}

.btn-action {
  margin: 0.25em;
}

.btn-dark {
  font-size: 15pt;
}

.btn-outline-dark {
  font-size: 15pt;
}

.alert-secondary {
  background-color: white;
  border-color: white;
  font-size: 1.4em;
  text-align: left;
}
</style>
