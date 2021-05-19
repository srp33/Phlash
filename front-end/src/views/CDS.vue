<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
      :logout="true"
    />
    <div class="container">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <h1>ID: {{ $route.params.cdsID }}</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <div>
          <div v-if="!viewOnly">Here you can make any necessary changes to this CDS by editing its
          start and stop sites and its function. Before making any final
          decisions you should review and consider all of the information below
          including the alternative open reading frames, the coding potential
          graphs, and the BLAST results. Click 'Save' to save changes made to this CDS including the notes 
          or you can click 'Delete' if you wish to remove this CDS.</div>
          <div v-else>Here you can view alternative open reading frames, the coding potential
          graphs, and the BLAST results. Click 'Save Notes' to save any edits that are made to the notes.
          The range of alternate open reading frames and BLAST results shown can be changed in 'Settings'.</div>
          Below you can see generic information for the current CDS including
          which auto-annotation programs called that gene.
        </div>
        <hr />
        <div style="float: right; width: 50%">
          <p><strong>Notes:</strong></p>
          <b-form-textarea
            id="textarea"
            v-model="notes"
            placeholder="Original Glimmer call @bp 408 has strength..."
            rows="5"
            max-rows="10"
          ></b-form-textarea>
        </div>
        <p>
          <strong>Left:</strong> {{ currentCDS.left }}<br />
          <strong>Right:</strong> {{ currentCDS.right }}<br />
          <span v-if="currentCDS.strand === '-'"
            ><strong>Strand:</strong> Direct</span
          >
          <span v-else><strong>Strand:</strong> Complementary</span><br />
          <strong>Frame:</strong> {{ this.frame }}<br />
          <strong>Called by:</strong> {{ this.calledBy }}<br />
          <strong>Product:</strong> {{ displayFunction }}
        </p>
        <button v-if="!viewOnly" type="button" class="btn btn-dark btn-action" @click="editCDS">
          <strong>&#9998; Save</strong>
        </button>
        <button v-else type="button" class="btn btn-dark btn-action" @click="editNotes">
          <strong>&#9998; Save Notes</strong>
        </button>
        <button
          v-if="!viewOnly"
          type="button"
          class="btn btn-dark btn-action"
          @click="deleteCDS($route.params.cdsID)"
        >
          <strong>&#128465; Delete</strong>
        </button>
        <hr />
        <div class="nav-btns-wrapper">
          <button
            type="button"
            class="btn btn-dark btn-nav"
            @click="navPrevCDS()"
          >
            <strong>&#129052; Prev CDS</strong>
          </button>
          <button
            type="button"
            class="btn btn-dark btn-nav"
            @click="navNextCDS()"
          >
            <strong>Next CDS &#129054;</strong>
          </button>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{
              name: 'Annotations',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129053; Return to Annotations</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
      <div class="coding-potential-table">
        <h4 style="text-align: center; margin: 1em">
          Alternative Open Reading Frames
        </h4>
        <div style="overflow: hidden">
          <div style="float: left; width: 40%">
            <strong style="font-size: 1.4em">Direct Strand</strong>
            <div class="table-responsive">
              <table id="cp-table" class="table table-hover">
                <tbody>
                  <tr v-for="(left, index) in dirLeftOptions" :key="index">
                    <th>{{ left }}-{{ dirRightOptions[index] }}</th>
                    <td>
                      <button
                        v-if="
                          left + dirRightOptions[index] !==
                          currentCDS.left + currentCDS.right
                        "
                        class="btn btn-dark btn-sm"
                        @click="setORF(left, dirRightOptions[index], '+')"
                      >
                        <strong v-if="!viewOnly">Select</strong>
                        <strong v-else>View</strong>
                      </button>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
          <div style="float: right; width: 40%">
            <strong style="font-size: 1.4em">Complementary Strand</strong>
            <div class="table-responsive">
              <table id="cp-table" class="table table-hover">
                <tbody>
                  <tr v-for="(left, index) in compLeftOptions" :key="index">
                    <th>{{ left }}-{{ compRightOptions[index] }}</th>
                    <td>
                      <button
                        v-if="
                          left + compRightOptions[index] !==
                          currentCDS.left + currentCDS.right
                        "
                        class="btn btn-dark btn-sm"
                        @click="setORF(left, compRightOptions[index], '-')"
                      >
                        <strong v-if="!viewOnly">Select</strong>
                        <strong v-else>View</strong>
                      </button>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
        </div>
      </div>
      <hr />
      <div class="coding-potential">
        <h4 style="text-align: center; margin: 1em; height: 100%">
          GeneMark's Coding Potential Per Frame
        </h4>
        <div class="alert alert-secondary">
          <p><strong>Key</strong></p>
          <p>
            <strong style="color: #2559aa">Blue Line:</strong> The coding
            potential in relation to base number.<br />
            <strong style="color: #d95f02">Orange Dashed Line:</strong> The
            selected start and stop positions.<br />
            <strong style="color: #1b9e77">Green Line:</strong> A 0.75 coding
            potential reference line.<br />
            <strong style="color: #7570b3">Purple Line:</strong> The previous
            gene's stop position and the next gene's start position.<br />
          </p>
        </div>
        <div class="coding-potential-graphs">
          <div v-if="dataExists">
            <Graphs
              :data1="data1"
              :data2="data2"
              :data3="data3"
              :data4="data4"
              :data5="data5"
              :data6="data6"
              :left="currentCDS.left"
              :right="currentCDS.right"
              :frame="frame"
              :prevRight="prevRight"
              :nextLeft="nextLeft"
            />
          </div>
          <div v-else>
            <p>Loading graphs...</p>
          </div>
        </div>
        <p class="graphs-caption">
          <em>Frames are counted 1-6 (direct 1-3 and complementary 4-6).</em>
        </p>
      </div>
      <hr />
      <h4 style="text-align: center; margin: 1em">
        BLAST Results for {{ newLeft }} - {{ newRight }}
      </h4>
      <BlastResults
        :blastResults="
          allBlastResults[
            currentCDS.left.toString() +
              '-' +
              currentCDS.right.toString() +
              '  ' +
              currentCDS.strand
          ]
        "
        :viewOnly="viewOnly"
        @newFunction="setFunction"
      />
      <hr />
      <div class="alert alert-secondary">
        <div style="float: right; width: 50%">
          <p><strong>Notes:</strong></p>
          <b-form-textarea
            id="textarea"
            v-model="notes"
            placeholder="Original Glimmer call @bp 408 has strength..."
            rows="5"
            max-rows="10"
          ></b-form-textarea>
        </div>
        <p>
          <strong>Left:</strong> {{ currentCDS.left }}<br />
          <strong>Right:</strong> {{ currentCDS.right }}<br />
          <span v-if="currentCDS.strand === '-'"
            ><strong>Strand:</strong> Direct</span
          >
          <span v-else><strong>Strand:</strong> Complementary</span><br />
          <strong>Frame:</strong> {{ this.frame }}<br />
          <strong>Called by:</strong> {{ this.calledBy }}<br />
          <strong>Product:</strong> {{ displayFunction }}
        </p>
        <button v-if="!viewOnly" type="button" class="btn btn-dark btn-action" @click="editCDS">
          <strong>&#9998; Save</strong>
        </button>
        <button v-else type="button" class="btn btn-dark btn-action" @click="editNotes">
          <strong>&#9998; Save Notes</strong>
        </button>
        <button
          v-if="!viewOnly"
          type="button"
          class="btn btn-dark btn-action"
          @click="deleteCDS($route.params.cdsID)"
        >
          <strong>&#128465; Delete</strong>
        </button>
        <hr />
        <div class="nav-btns-wrapper">
          <button
            type="button"
            class="btn btn-dark btn-action"
            @click="navPrevCDS()"
          >
            <strong>&#129052; Prev CDS</strong>
          </button>
          <button
            type="button"
            class="btn btn-dark btn-action"
            @click="navNextCDS()"
          >
            <strong>Next CDS &#129054;</strong>
          </button>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{
              name: 'Annotations',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129053; Return to Annotations</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
    <b-toast id="cds-status" variant="primary">
      <template #toast-title>
        <strong class="text-size"> {{statusTitle}} </strong>
      </template>
      <div class="text-size">{{ statusMessage }}</div>
    </b-toast>
  </div>
</template>

<script>
import axios from 'axios';
import Navbar from '../components/Navbar.vue';
import BlastResults from '../components/BlastResults.vue';
import Graphs from '../components/Graphs.vue';
import Loading from 'vue-loading-overlay';
import 'vue-loading-overlay/dist/vue-loading.css';
import Vue from 'vue';

export default {
  name: 'CDS',
  components: {
    BlastResults,
    Graphs,
    Loading,
    Navbar,
  },

  data() {
    return {
      viewOnly: false,
      dirLeftOptions: [],
      compLeftOptions: [],
      dirRightOptions: [],
      compRightOptions: [],
      dirBlastResults: [],
      compBlastResults: [],
      allBlastResults: [],
      currentCDS: {
        id: '',
        left: '',
        right: '',
        strand: '',
        function: '',
        status: '',
        frame: '',
      },
      updatedCDS: {
        id: '',
        left: '',
        right: '',
        strand: '',
        function: '',
        status: '',
      },
      frame: null,
      newFunction: 'None selected',
      displayFunction: '',
      newLeft: null,
      newRight: null,
      newStrand: null,
      data1: [{ x: [], y: [] }],
      data2: [{ x: [], y: [] }],
      data3: [{ x: [], y: [] }],
      data4: [{ x: [], y: [] }],
      data5: [{ x: [], y: [] }],
      data6: [{ x: [], y: [] }],
      nextCDS: null,
      prevCDS: null,
      nextLeft: null,
      prevRight: null,
      calledBy: '',
      glimmer: '',
      genemark: '',
      phanotate: '',
      notes: '',
      dataExists: false,
      pageLoading: true,
      showFunction: false,
      showLeft: false,
      saved: true,
      statusMessage: "",
      statusTitle: "",
    };
  },

  beforeCreate() {
    Vue.GoogleAuth.then(auth2 => {
      if (!auth2.isSignedIn.get()) {
        this.$router.push('/');
      }
      axios
        .get(process.env.VUE_APP_BASE_URL + `/check_user/${auth2.currentUser.get().getBasicProfile().getEmail()}/${this.$route.params.phageID}`)
        .then((response) => {
          if (response.data === "fail") {
            this.$router.push('/');
          }
          else if (response.data.view) {
            this.viewOnly = true
          }
        })
        .catch((error) => {
          console.error(error);
        });
    })
  },

  created() {
    this.getData(this.$route.params.cdsID);
    window.scrollTo(0, 0);
  },

  computed: {
    dirBlastKeys: function () {
      return Object.keys(this.dirBlastResults);
    },

    compBlastKeys: function () {
      return Object.keys(this.compBlastResults);
    },

    navUpload: function () {
      return !this.viewOnly;
    },

    navBlast: function () {
      return !this.viewOnly;
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
     * Gets the GeneMark, Blast, and auto-annotation data for the given CDS.
     * @param {string} cdsID the ID of the CDS.
     */
    getData(cdsID) {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${cdsID}`
        )
        .then((response) => {
          if (response.data.message !== 'Finished') {
            this.$router.push(`/annotations/${this.$route.params.phageID}`);
          }
          this.currentCDS = response.data.cds;
          this.updatedCDS.id = this.currentCDS.id;
          this.updatedCDS.left = this.currentCDS.left;
          this.updatedCDS.right = this.currentCDS.right;
          this.updatedCDS.strand = this.currentCDS.strand;
          this.updatedCDS.status = this.currentCDS.status;
          this.updatedCDS.function = this.currentCDS.function;
          if (this.currentCDS.function !== '@DELETED') {
            this.newFunction = this.currentCDS.function;
            this.displayFunction = this.newFunction;
            let indexSeparation = this.newFunction.indexOf('##');
            if (indexSeparation !== -1) {
              this.displayFunction = this.displayFunction.substring(
                0,
                indexSeparation
              );
            }
            if (this.newFunction[0] === '@') {
              this.displayFunction = this.displayFunction.substring(1);
            }
          }
          this.dirBlastResults = response.data.dir_blast;
          this.compBlastResults = response.data.comp_blast;
          this.allBlastResults = response.data.all_blast;
          this.dirLeftOptions = response.data.dir_left_options;
          this.dirRightOptions = response.data.dir_right_options;
          this.compLeftOptions = response.data.comp_left_options;
          this.compRightOptions = response.data.comp_right_options;
          this.newLeft = this.currentCDS.left;
          this.newRight = this.currentCDS.right;
          this.newStrand = this.currentCDS.strand;
          this.prevRight = response.data.prev_right;
          this.nextLeft = response.data.next_left;
          this.data1 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_1,
            },
          ];
          this.data2 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_2,
            },
          ];
          this.data3 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_3,
            },
          ];
          this.data4 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_4,
            },
          ];
          this.data5 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_5,
            },
          ];
          this.data6 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_6,
            },
          ];
          this.dataExists = true;
          this.pageLoading = false;
          this.frame = this.currentCDS.frame;
          this.notes = this.currentCDS.notes;
          this.glimmer = response.data.glimmer;
          this.genemark = response.data.genemark;
          this.phanotate = response.data.phanotate;
          var called = false;
          var cds =
            this.currentCDS.left.toString() +
            '-' +
            this.currentCDS.right.toString() +
            ' ' +
            this.currentCDS.strand;
          if (this.glimmer.indexOf(cds) > -1) {
            called = true;
            this.calledBy += 'Glimmer, ';
          }
          if (this.genemark.indexOf(cds) > -1) {
            called = true;
            this.calledBy += 'GeneMark, ';
          }
          if (this.phanotate.indexOf(cds) > -1) {
            called = true;
            this.calledBy += 'Phanotate, ';
          }
          if (called) {
            this.calledBy = this.calledBy.substring(
              0,
              this.calledBy.length - 2
            );
          } else {
            this.calledBy = 'None';
          }
          this.nextCDS = response.data.nextCDS;
          this.prevCDS = response.data.prevCDS;
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Re-routes to the next CDS and/or displays a save reminder.
     */
    navNextCDS() {
      var cont = true;
      if (!this.saved && !this.viewOnly) {
        cont = confirm('Are you sure you want to continue without saving?');
      }
      if (cont === true) {
        if (this.nextCDS !== 'undefined') {
          this.$route.params.cdsID = this.nextCDS;
          this.$router.push(
            `/annotations/cds/${this.$route.params.phageID}/${this.$route.params.cdsID}`
          );
          window.location.reload();
        } else {
          this.$router.push(`/annotations/${this.$route.params.phageID}`);
        }
      }
    },

    /**
     * Re-routes to the previous CDS and/or displays a save reminder.
     */
    navPrevCDS() {
      var cont = true;
      if (!this.saved && !this.viewOnly) {
        cont = confirm('Are you sure you want to continue without saving?');
      }
      if (cont === true) {
        if (this.prevCDS !== 'undefined') {
          this.$route.params.cdsID = this.prevCDS;
          this.$router.push(
            `/annotations/cds/${this.$route.params.phageID}/${this.$route.params.cdsID}`
          );
          window.location.reload();
        } else {
          this.$router.push(`/annotations/${this.$route.params.phageID}`);
        }
      }
    },

    /**
     * Changes the CDS data to reflect the user's changes.
     */
    editCDS() {
      this.saved = true;
      this.updatedCDS = this.currentCDS;
      this.updatedCDS.left = this.newLeft;
      this.updatedCDS.right = this.newRight;
      if (this.newFunction === 'None selected') {
        this.newFunction = '@' + this.newFunction;
      }
      const payload = {
        id: this.updatedCDS.id,
        left: this.updatedCDS.left,
        right: this.updatedCDS.right,
        strand: this.updatedCDS.strand,
        function: this.newFunction,
        notes: this.notes,
        status: this.currentCDS.status,
      };
      this.updateCDS(payload, this.updatedCDS.id);
    },

    /**
     * Changes the CDS data to reflect the user's changes.
     */
    editNotes() {
      const payload = {
        notes: this.notes,
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${this.updatedCDS.id}`,
          payload
        )
        .then(() => {
          this.statusMessage = `The CDS ${this.updatedCDS.id} has been saved.`;
          this.statusTitle = "SAVED";
          this.$bvToast.show('cds-status');
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Updates the database with the new changes on the current cds.
     * @param {dictionary} payload the new CDS data.
     * @param {string} cdsID the ID of the CDS to be updated.
     */
    updateCDS(payload, cdsID) {
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${cdsID}`,
          payload
        )
        .then(() => {
          if (this.newFunction !== '@DELETED') {
            this.statusMessage = `The CDS ${cdsID} has been saved.`;
            this.statusTitle = "SAVED";
            this.$bvToast.show('cds-status');
          } else {
            this.statusMessage = `The CDS ${cdsID} has been deleted. You will be re-routed to the next CDS.`;
            this.statusTitle = "DELETED";
            this.$bvToast.show('cds-status')
            this.navNextCDS();
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Changes the function of the CDS to deleted.
     * @param {string} cdsID the ID of the CDS to be deleted.
     */
    deleteCDS(cdsID) {
      this.newFunction = '@DELETED';
      this.editCDS();
    },

    /**
     * Updates the function.
     * @param {string} function the user selected function.
     */
    setFunction(funct) {
      this.saved = false;
      let match = funct.match(/(.*)##(.*)/);
      this.displayFunction = match[1];
      this.newFunction = '@' + funct;
    },

    /**
     * Updates the left and right postions to what the user selected.
     * @param {number} left the user selected left.
     * @param {number} right the user selected right.
     */
    setORF(left, right, strand) {
      this.saved = false;
      if (left !== this.newLeft || right !== this.newRight) {
        this.dataExists = false;
        this.newFunction = '@None selected';
        this.displayFunction = 'None selected';
        this.newRight = right;
        this.currentCDS.right = right;
        this.newLeft = left;
        this.newStrand = strand;
        this.currentCDS.left = left;
        this.currentCDS.strand = strand;
        this.calledBy = '';
        var called = false;
        var cds =
          this.currentCDS.left.toString() +
          '-' +
          this.currentCDS.right.toString() +
          ' ' +
          this.currentCDS.strand;
        if (this.glimmer.indexOf(cds) > -1) {
          called = true;
          this.calledBy += 'Glimmer, ';
        }
        if (this.genemark.indexOf(cds) > -1) {
          called = true;
          this.calledBy += 'GeneMark, ';
        }
        if (this.phanotate.indexOf(cds) > -1) {
          called = true;
          this.calledBy += 'Phanotate, ';
        }
        if (called) {
          this.calledBy = this.calledBy.substring(0, this.calledBy.length - 2);
        } else {
          this.calledBy = 'None';
        }
        this.frame = ((left + 2) % 3) + 1;
        if (this.currentCDS.strand === '-') {
          this.frame = ((this.currentCDS.right + 2) % 3) + 4;
        } else this.frame = ((this.currentCDS.left + 2) % 3) + 1;
        this.$nextTick().then(() => {
          this.dataExists = true;
        });
      }
    },
  },
};
</script>

<style scoped>
.blast-results {
  margin: 1em;
  margin-bottom: 3em;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 0.25em;
}

h1 {
  margin-top: 0.7em;
}

h4 {
  font-size: 1.7em;
}

.info-bottom {
  margin: 1em auto;
}

.btn-action {
  margin: 0.25em;
}

.coding-potential-table {
  margin: 2em auto;
  width: 80%;
}

.table-responsive {
  max-height: 15em;
  overflow-y: auto;
  display: inline-block;
  font-size: 1.4em;
}

.table-responsive thead th {
  position: sticky;
  top: 0;
  background: #eee;
  border: darkgray;
  font-size: 1em;
}

caption {
  display: table-caption;
}

tbody {
  width: 100%;
}

.coding-potential-graphs {
  height: 40em;
}

.graphs-caption {
  margin-top: 1em;
  text-align: center;
  color: grey;
}

.alert-secondary {
  background-color: white;
  border-color: white;
  font-size: 1.4em;
  text-align: left;
}

.btn-dark {
  font-size: 15pt;
}

.text-size {
  font-size: 1.2em;
}

@media only screen and (max-width: 50em) {
  .coding-potential-graphs {
    height: 50em;
  }
}
</style>
