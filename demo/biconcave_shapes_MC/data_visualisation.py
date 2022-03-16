import json
from typing import Union, List

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import product
from matplotlib.colors import LogNorm


def _normed_vec(vec: np.ndarray):
    norm = vec.dot(vec)
    if norm > 0:
        return vec/np.sqrt(norm)
    else:
        return np.zeros(vec.size)


def _face_name_splitter(face_name: str):
    return np.array(list(map(int, face_name.split('_'))))


def _edge_namer(id1: int, id2: int) -> str:
    if id1 < id2:
        edge_list = [id1, id2]
    else:
        edge_list = [id2, id1]
    return ''.join(map(str, edge_list))


def find_two_neighbours(node_nns: List[int], nn_nns: List[int]):
    return list(set(node_nns).intersection(nn_nns))


def triangle_id_maker(node_id:int, nn_id:int, cnn_id:int):
    ids = [node_id, nn_id, cnn_id]
    ids.sort()
    return '_'.join(list(map(str,ids)))


class SingleTimeFrameImporter(object):

    def __init__(self, file_name: str):
        with open(file_name, "r") as f:
            self.data_dict: dict = json.load(f)

    @classmethod
    def load_data(cls, file_name: str):
        instance = cls(file_name)
        instance._load_egg()
        instance._load_energy_data()
        instance._make_geometry()
        return instance

    @classmethod
    def egg_import(cls, file_name: str):
        instance =cls(file_name)
        instance._load_egg()
        instance._make_geometry()
        return instance

    def _load_egg(self):
        node_dict = self.data_dict
        self.node_dict = dict(sorted(zip(list(map(int, node_dict.keys())), node_dict.values())))
        self.no_cell = True

    def _load_energy_data(self):
        pass

    def _make_geometry(self):
        self.overlapers = []
        self._make_positions()
        self._make_edges()
        self._build_triangles()
        self._make_faces()

    def xyz_export_data(self) -> str:

        data_array = np.array([np.array(self.xs, dtype=str),
                               np.array(self.ys, dtype=str),
                               np.array(self.zs, dtype=str)]).T
        if self.no_cell:
            n = len(self.xs)
            text: str = f'{n}\nAtoms {self._out_stream_constructor("1", data_array)}'
        else:
            cell_data_array = np.array([np.array(self.cell_xs, dtype=str),
                                        np.array(self.cell_ys, dtype=str),
                                        np.array(self.cell_zs, dtype=str)]).T
            n = len(self.xs) + len(self.cell_xs)
            text: str = f'{n}\nAtoms {self._out_stream_constructor("1", data_array)}'
        return text

    @staticmethod
    def _out_stream_constructor(name_stream: Union[str, List[str]], xyz_matrix):
        text = ''
        if type(name_stream) is str:
            name_stream = [name_stream for _ in range(len(xyz_matrix))]
        for name, row in zip(name_stream, xyz_matrix):
            line: str = f"{name} {' '.join(list(map(str, row)))}"
            text = '\n'.join([text, line])
        return text

    def xyz_save_to_file(self, file_name):
        with open(file_name, 'a') as f:
            f.write(self.xyz_export_data())
        return f

    def xyz_append_to_file(self, f):
        f.write(self.xyz_export_data())
        return f

    def exyz_save_to_file(self, file_name):
        with open(file_name, 'a') as f:
            f.write(self.exyz_export_data())
        return f

    def exyz_append_to_file(self, f):
        f.write(self.exyz_export_data())
        return f

    def _make_positions(self):
        xs = []
        ys = []
        zs = []
        cvs = []

        for (k, v) in sorted(self.node_dict.items()):
            xs.append(v['pos'][0])
            ys.append(v['pos'][1])
            zs.append(v['pos'][2])
            curv = np.array(v['curvature_vec'])
            norm = np.sqrt(curv.dot(curv))
            cvs.append(norm)

        self.xs = np.array(xs)
        self.ys = np.array(ys)
        self.zs = np.array(zs)
        self.node_curvatures = np.array(cvs)

        if not self.no_cell:
            for (k, v) in sorted(self.cell_node_dict.items()):
                xs.append(v['pos'][0])
                ys.append(v['pos'][1])
                zs.append(v['pos'][2])
                curv = np.array(v['curvature_vec'])
                norm = np.sqrt(curv.dot(curv))
                cvs.append(norm)

            self.cell_xs = np.array(xs)
            self.cell_ys = np.array(ys)
            self.cell_zs = np.array(zs)
            self.cell_node_curvatures = np.array(cvs)

    def _make_edges(self):
        self.edge_dict = {}
        lines = []
        for node_id, node in self.node_dict.items():
            for nn_id in node['nn_ids']:
                edge_id = _edge_namer(node_id, nn_id)
                if edge_id not in self.edge_dict:
                    line = np.array([self.node_dict[node_id]['pos'], self.node_dict[nn_id]['pos']])
                    lines.append(line)
                    self.edge_dict[edge_id] = {'nodes': [node_id, nn_id], 'line': line}
        self.lines = np.array(lines)

    def _build_triangles(self):
        for node_id, node in self.node_dict.items():
            node["triangle_ids"] = set()
            nn_ids = node["nn_ids"]
            for nn_id in nn_ids:
                two_neighbours = find_two_neighbours(node["nn_ids"], self.node_dict[nn_id]["nn_ids"])
                assert len(two_neighbours)==2
                node["triangle_ids"].add(triangle_id_maker(node_id, nn_id, two_neighbours[0]))
                node["triangle_ids"].add(triangle_id_maker(node_id, nn_id, two_neighbours[1]))

            node["triangle_ids"] = list(node["triangle_ids"])

    def _make_faces(self):
        faces = []
        face_curves = []
        face_names = {}
        face_curve_max = 0
        face_curve_min = 1e5

        for _faces in self.node_dict.values():
            for _face in _faces['triangle_ids']:
                if _face not in face_names:
                    face_names[_face] = 1
                    n0, n1, n2 = _face_name_splitter(_face)
                    p0 = np.array(self.node_dict[n0]['pos'])
                    p1 = np.array(self.node_dict[n1]['pos'])
                    p2 = np.array(self.node_dict[n2]['pos'])

                    faces.append([p0, p1, p2])

                    # try:
                    k0 = np.array(self.node_dict[n0]['curvature_vec'])
                    k1 = np.array(self.node_dict[n1]['curvature_vec'])
                    k2 = np.array(self.node_dict[n2]['curvature_vec'])
                    face_norm = (k0+k1+k2)/3
                    face_curve = np.sqrt(face_norm.dot(face_norm))

                    if face_curve < face_curve_min:
                        face_curve_min=face_curve
                    if face_curve > face_curve_max:
                        face_curve_max = face_curve
                    face_curves.append(face_curve)

        self.faces = np.array(faces)
        self.face_curve_max = face_curve_max
        self.face_curve_min = face_curve_min
        if face_curve_max!=0:
            self.face_curves = np.array(face_curves)/face_curve_max
        else:
            self.face_curves = np.array(face_curves)

    def build_polygon(self, alpha):
        cmap = [cm.get_cmap('Reds')(c) for c in self.face_curves]
        return Poly3DCollection(self.faces,
                         edgecolors='k', linewidths=0.5, facecolors=cmap, alpha=alpha,
                         norm=LogNorm(vmin=self.face_curve_min, vmax=self.face_curve_max)
                         ), cmap


    def make_polygon(self, stretch=1.1, **kwargs):
        no_box = kwargs.get("no_box", True)
        alpha = kwargs.get("alpha", 1)
        self._make_empty_figure()
        polygons, cmap = self.build_polygon(alpha=alpha)
        self.ax.add_collection3d(polygons)

        if no_box:
            self.ax.set_axis_off()

        cmp = np.round(np.linspace(self.face_curve_min, self.face_curve_max, 20), 3)
        z = cmp[:, None]
        self.cax.tick_params(axis='y', labelright=True, right=True, labelleft=False, left=False)
        self.cax.imshow(z, cmap="Reds", origin='lower')
        self.cax.set_label('local mean curv.')
        self.cax.set_yticks(np.arange(0, cmp.size))
        self.cax.set_yticklabels([str(tick) for tick in cmp])
        self.cax.set_xticks([])
        all_min = np.min([self.xs, self.ys, self.zs])
        all_max = np.max([self.xs, self.ys, self.zs])
        self.ax.set_xlim(all_min, all_max)
        self.ax.set_ylim(all_min, all_max)
        self.ax.set_zlim(all_min, all_max)

        return self


    def _make_empty_figure(self):
        spec = gridspec.GridSpec(ncols=2, nrows=1,
                                 width_ratios=[10, 1], wspace=0.1,
                                 hspace=0.1)
        self.fig: plt.Figure = plt.figure()
        self.ax: plt.Axes = self.fig.add_subplot(spec[0], projection='3d')
        self.cax: plt.Axes = self.fig.add_subplot(spec[1])

    def save_plot(self, save_path: str, **kwargs):
        self.fig.savefig(save_path, **kwargs)

    def show_plot(self):
        self.fig.show()


class TopologyChecks(object):
    def __init__(self, single_tf_data_set: SingleTimeFrameImporter, n_it):
        self.n_it = n_it
        self.local_dataset: SingleTimeFrameImporter = single_tf_data_set
        self.n_nodes = len(self.local_dataset.node_dict)
        self.n_edges = len(self.local_dataset.lines)
        self.n_edges2 = sum([len(self.local_dataset.node_dict[i]["nn_ids"]) for i in range(self.n_nodes)])/2
        self.n_faces = len(self.local_dataset.faces)
        self.n_nodes_exp = 12 + 30 * (2 ** self.n_it - 1) + 20 * (2 ** (self.n_it - 1) * (2 ** self.n_it - 3) + 1)
        self.n_edges_exp = 30 * (4 ** self.n_it)
        self.n_faces_exp = 20 * (4 ** self.n_it)

    def check_euler_characteristic(self):
        print("N_NODES: ", self.n_nodes, "N_EDGES:", self.n_edges, "or:", self.n_edges2, "N_FACES: ", self.n_faces)
        print("expected")
        print("N_NODES: ", self.n_nodes_exp, "N_EDGES: ", self.n_edges_exp, "N_FACES: ", self.n_faces_exp)
        print("Euler Characteristic: ", self.n_nodes - self.n_edges + self.n_faces)

    def check_node_name_contiguity(self):
        k_max = max(self.local_dataset.node_dict.keys())
        print("N_IT:", self.n_it, "\n Max_id + 1 - N_NODES: \n",
              "counter mismatch: ", k_max + 1 - self.n_nodes, '\n')

    def check_all(self):
        print("checking")
        self.check_euler_characteristic()
        self.check_node_name_contiguity()
        print(80*'=')


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    import json

    egg_init = SingleTimeFrameImporter.egg_import(f'test_run_init.json')
    egg_init.make_polygon()
    egg_final = SingleTimeFrameImporter.egg_import(f'test_run_final.json')
    egg_final.make_polygon(alpha=0.8)
    plt.show()
