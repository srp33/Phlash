import Vue from 'vue';
import Router from 'vue-router';
import Home from './views/Home.vue';
import Upload from './views/Upload.vue';
import DNAMaster from './views/DNAMaster.vue';
import Blast from './views/Blast.vue';
import Annotations from './views/Annotations.vue';
import Fail from './views/Fail.vue';
import More from './views/More.vue';
import Pass from './views/Pass.vue';
import GenBank from './views/GenBank.vue';


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
         path: '/upload/:currentUser',
         name: 'Upload',
         component: Upload,
      },
      {
         path: '/dnamaster/:currentUser',
         name: 'DNAMaster',
         component: DNAMaster,
      },
      {
         path: '/blast/:currentUser',
         name: 'Blast',
         component: Blast
      },
      {
         path: '/annotations/:currentUser',
         name: 'Annotations',
         component: Annotations,
      },
      {
         path: '/annotations/pass/:currentUser/:id',
         name: 'Pass',
         component: Pass,
      },
      {
         path: '/annotations/fail/:currentUser/:id',
         name: 'Fail',
         component: Fail,
      },
      {
         path: '/annotations/more/:currentUser/:id',
         name: 'More',
         component: More,
      },
      {
         path: '/genbank',
         name: 'GenBank',
         component: GenBank,
      },
   ],
});
