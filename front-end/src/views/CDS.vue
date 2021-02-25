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
        <h1>ID: {{ $route.params.cdsID }}</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p>
          Here you can make any necessary changes to this CDS by editing its start and stop sites and its function. 
          Before making any final decisions you should review and consider all of the information below including 
          the alternative open reading frames, the coding potential graphs, and the BLAST results. 
          Click 'Save' to save changes made to this CDS or you can click 'Delete' if you wish to remove this CDS. 
          Above you can see generic information for the current CDS including which auto-annotation programs called that gene.
        </p>
        <hr />
        <div style="float: right; width: 50%;">
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
          <strong>Left:</strong> {{ currentCDS.start }}<br />
          <strong>Right:</strong> {{ currentCDS.stop }}<br />
          <span v-if="currentCDS.strand == '-'"><strong>Strand:</strong> Direct</span>
          <span v-else><strong>Strand:</strong> Complementary</span><br />
          <strong>Frame:</strong> {{ this.frame }}<br />
          <strong>Called by:</strong> {{ this.calledBy }}<br />
          <strong>Product:</strong> {{displayFunction}}
        </p>
        <button
          type="button"
          class="btn btn-dark btn-action"
          @click="editCDS"
        >
          <strong>&#9998; Save</strong>
        </button>
        <button
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
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129053; Return to Annotations</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
      <div class="coding-potential-table">
        <h4 style="text-align: center; margin: 1em">Alternative Open Reading Frames</h4>
        <div style="overflow: hidden;">
          <div style="float: left; width: 40%;">
            <strong style="font-size:1.4em;">Direct Strand</strong>
            <div class="table-responsive">
              <table id="cp-table" class="table table-hover">
                <tbody>
                  <tr v-for="(start, index) in dirStartOptions" :key="index">
                    <th>{{ start }}-{{ dirStopOptions[index] }}</th>
                    <td>
                      <button
                        v-if="start + dirStopOptions[index] !== currentCDS.start + currentCDS.stop"
                        class="btn btn-dark btn-sm"
                        @click="setORF(start, dirStopOptions[index], '+')"
                      >
                        <strong>Select</strong>
                      </button>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
          <div style="float: right; width: 40%">
            <strong style="font-size:1.4em;">Complementary Strand</strong>
            <div class="table-responsive">
              <table id="cp-table" class="table table-hover">
                <tbody>
                  <tr v-for="(start, index) in compStartOptions" :key="index">
                    <th>{{ start }}-{{ compStopOptions[index] }}</th>
                    <td>
                      <button
                        v-if="start + compStopOptions[index] !== currentCDS.start + currentCDS.stop"
                        class="btn btn-dark btn-sm"
                        @click="setORF(start, compStopOptions[index], '-')"
                      >
                        <strong>Select</strong>
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
            <strong style="color:#2559AA;">Blue Line:</strong> The coding potential in relation to base number.<br />
            <strong style="color:#d95f02;">Orange Dashed Line:</strong> The selected start and stop positions.<br />
            <strong style="color:#1b9e77;">Green Line:</strong> A 0.75 coding potential reference line.<br />
            <strong style="color:#7570b3;">Purple Line:</strong> The previous gene's stop position and the next
            gene's start position.<br />
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
              :start="currentCDS.start"
              :stop="currentCDS.stop"
              :frame="frame"
              :prevStop="prevStop"
              :nextStart="nextStart"
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
      <h4 style="text-align: center; margin: 1em">BLAST Results for {{ newStart }} - {{ newStop }}</h4>
      <BlastResults
        :blastResults="allBlastResults[currentCDS.start.toString() + '-' + currentCDS.stop.toString() + '  ' + currentCDS.strand]"
        :allowSelect="true"
        @newFunction="setFunction"
      />
      <hr />
      <div class="alert alert-secondary">
        <div style="float: right; width: 50%;">
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
          <strong>Left:</strong> {{ currentCDS.start }}<br />
          <strong>Right:</strong> {{ currentCDS.stop }}<br />
          <span v-if="currentCDS.strand == '-'"><strong>Strand:</strong> Direct</span>
          <span v-else><strong>Strand:</strong> Complementary</span><br />
          <strong>Frame:</strong> {{ this.frame }}<br />
          <strong>Called by:</strong> {{ this.calledBy }}<br />
          <strong>Product:</strong> {{displayFunction}}
        </p>
        <button
          type="button"
          class="btn btn-dark btn-action"
          @click="editCDS"
        >
          <strong>&#9998; Save</strong>
        </button>
        <button
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
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129053; Return to Annotations</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue";
import BlastResults from "../components/BlastResults.vue";
import Graphs from "../components/Graphs.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "CDS",
  components: {
    BlastResults,
    Graphs,
    Loading,
    Navbar,
  },

  data() {

    return {
      dirStartOptions: [],
      compStartOptions: [],
      dirStopOptions: [],
      compStopOptions: [],
      dirBlastResults: [],
      compBlastResults: [],
      allBlastResults: [],
      currentCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: "",
        frame: "",
      },
      updatedCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: "",
      },
      frame: null,
      newFunction: "None selected",
      displayFunction: "",
      newStart: null,
      newStop: null,
      newStrand: null,
      data1: [{ x: [], y: [] }],
      data2: [{ x: [], y: [] }],
      data3: [{ x: [], y: [] }],
      data4: [{ x: [], y: [] }],
      data5: [{ x: [], y: [] }],
      data6: [{ x: [], y: [] }],
      nextCDS: null,
      prevCDS: null,
      nextStart: null,
      prevStop: null,
      calledBy: "",
      glimmer: "",
      genemark: "",
      phanotate: "",
      notes: "",
      dataExists: false,
      pageLoading: true,
      showFunction: false,
      showStart: false,
      saved: true,
    };

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
          if (response.data.message != "Finished") {
            this.$router.push(
              `/annotations/${this.$route.params.phageID}`
            );
          }
          this.currentCDS = response.data.cds;
          if (this.currentCDS.function != "DELETED") {
            this.newFunction = this.currentCDS.function;
            this.displayFunction = this.newFunction;
            let indexSeparation = this.newFunction.indexOf('##');
            if (indexSeparation != -1) {
              this.displayFunction = this.displayFunction.substring(0, indexSeparation);
            }
            if (this.newFunction[0] == '@') {
              this.displayFunction = this.displayFunction.substring(1);
            }
            console.log(this.displayFunction);
          }
          this.dirBlastResults = response.data.dir_blast;
          this.compBlastResults = response.data.comp_blast;
          this.allBlastResults = response.data.all_blast;
          this.dirStartOptions = response.data.dir_start_options;
          this.dirStopOptions = response.data.dir_stop_options;
          this.compStartOptions = response.data.comp_start_options;
          this.compStopOptions = response.data.comp_stop_options;
          this.newStart = this.currentCDS.start;
          this.newStop = this.currentCDS.stop;
          this.newStrand = this.currentCDS.strand;
          this.prevStop = response.data.prev_stop;
          this.nextStart = response.data.next_start;
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
          var cds = this.currentCDS.start.toString() + '-' + this.currentCDS.stop.toString() + ' ' + this.currentCDS.strand;
          if (this.glimmer.indexOf(cds) > -1) {
            called = true;
            this.calledBy += "Glimmer, ";
          }
          if (this.genemark.indexOf(cds) > -1) {
            called = true;
            this.calledBy += "GeneMark, ";
          }
          if (this.phanotate.indexOf(cds) > -1) {
            called = true;
            this.calledBy += "Phanotate, ";
          }
          if (called) { this.calledBy = this.calledBy.substring(0, this.calledBy.length - 2); }
          else { this.calledBy = "None"}
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
      if (!this.saved) {
        cont = confirm("Are you sure you want to continue without saving?");
      }
      if (cont == true) {
        if (this.nextCDS != "undefined") {
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
      if (!this.saved) {
        cont = confirm("Are you sure you want to continue without saving?");
      }
      if (cont == true) {
        if (this.prevCDS != "undefined") {
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
      this.updatedCDS.start = this.newStart;
      this.updatedCDS.stop = this.newStop;
      this.updatedCDS.function = "@" + this.newFunction;
      const payload = {
        id: this.updatedCDS.id,
        start: this.updatedCDS.start,
        stop: this.updatedCDS.stop,
        strand: this.updatedCDS.strand,
        function: this.updatedCDS.function,
        notes: this.notes,
        status: this.currentCDS.status,
      };
      this.updateCDS(payload, this.updatedCDS.id);
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
          if (this.newFunction != "DELETED") {
            this.$bvToast.toast(
                `The CDS ${cdsID} has been saved.`,
                {
                  title: "SAVED",
                  autoHideDelay: 5000,
                  appendToast: false,
                }
              );
          }
          else {
            this.$bvToast.toast(
                `The CDS ${cdsID} has been deleted. You will be re-routed to the next CDS.`,
                {
                  title: "DELETED",
                  autoHideDelay: 5000,
                  appendToast: false,
                }
              );
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
      this.newFunction = "DELETED";
      this.editCDS();
    },

    /**
     * Updates the function.
     * @param {string} function the user selected function.
     */
    setFunction(funct) {
      console.log(funct);
      this.saved = false;
      let match = funct.match(/(.*)##(.*)/);
      this.displayFunction = match[1];
      console.log(this.displayFunction);
      this.newFunction = funct;
    },

    /**
     * Updates the start and stop postions to what the user selected.
     * @param {number} start the user selected start.
     * @param {number} stop the user selected stop.
     */
    setORF(start, stop, strand) {
      this.saved = false;
      if (start != this.newStart || stop != this.newStop) {
        this.dataExists = false;
        this.newFunction = "";
        this.displayFunction = "None selected";
        this.newStop = stop;
        this.currentCDS.stop = stop;
        this.newStart = start;
        this.newStrand = strand;
        this.currentCDS.start = start;
        this.currentCDS.strand = strand;
        this.calledBy = "";
        var called = false;
        var cds = this.currentCDS.start.toString() + '-' + this.currentCDS.stop.toString() + ' ' + this.currentCDS.strand;
        if (this.glimmer.indexOf(cds) > -1) {
          called = true;
          this.calledBy += "Glimmer, ";
        }
        if (this.genemark.indexOf(cds) > -1) {
          called = true;
          this.calledBy += "GeneMark, ";
        }
        if (this.phanotate.indexOf(cds) > -1) {
          called = true;
          this.calledBy += "Phanotate, ";
        }
        if (called) { this.calledBy = this.calledBy.substring(0, this.calledBy.length - 2); }
        else { this.calledBy = "None"}
        this.frame = ((start + 2) % 3) + 1;
        if (this.currentCDS.strand == "-") {
          this.frame = ((this.currentCDS.stop + 2) % 3) + 4;
        }
        else this.frame = ((this.currentCDS.start + 2) % 3) + 1;
        console.log(this.frame);
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
  margin-top: .7em;
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
  text-align: left;
  color: grey;
}

.alert-secondary {
  background-color: white;
  border-color:white;
  font-size: 1.40em;
  text-align: left;
}

.btn-dark {
  font-size: 15pt;
}

@media only screen and (max-width: 50em) {
  .coding-potential-graphs {
    height: 50em;
  }
}
</style>
