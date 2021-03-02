import Vue from 'vue';
import Router from 'vue-router';
import Home from './views/Home.vue';
import Upload from './views/Upload.vue';
import Blast from './views/Blast.vue';
import Annotations from './views/Annotations.vue';
import CDS from './views/CDS.vue';
import GeneMap from './views/GeneMap.vue';
import GenBank from './views/GenBank.vue';
import Contact from './views/Contact.vue';

Vue.use(Router);

export default new Router({
  mode: 'history',
  base: process.env.BASE_URL,
  routes: [
    {
      path: '/',
      name: 'Home',
      component: Home,
    },
    {
      path: '/upload/:phageID',
      name: 'Upload',
      component: Upload,
    },
    {
      path: '/blast/:phageID',
      name: 'Blast',
      component: Blast,
    },
    {
      path: '/annotations/:phageID',
      name: 'Annotations',
      component: Annotations,
    },
    {
      path: '/annotations/cds/:phageID/:cdsID',
      name: 'CDS',
      component: CDS,
    },
    {
      path: '/geneMap/:phageID',
      name: 'GeneMap',
      component: GeneMap,
    },
    {
      path: '/genBank/:phageID',
      name: 'GenBank',
      component: GenBank,
    },
    {
      path: '/contact',
      name: 'Contact',
      component: Contact,
    },
  ],
});
