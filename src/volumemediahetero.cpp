/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <nori/volumemedia.h>
#include <nori/volumeemitterhetero.h>

NORI_NAMESPACE_BEGIN

class VolumeMediaHetero : public VolumeMedia {
public:
    VolumeMediaHetero(const PropertyList &props) {
        m_name = props.getString("name", "");
        m_bsdf = (BSDF *)(NoriObjectFactory::createInstance("volumeHG", props));
        theta_t = props.getFloat("theta_t", 0.1f);
        Vector3f scale = props.getVector("scale", Vector3f(1.0f, 1.0f, 1.0f));

        Transform local_transform = props.getTransform("local", Transform());

        toWorld = props.getTransform("toWorld", Transform());
        toWorld = toWorld * local_transform;
        toLocal = toWorld.inverse();
        density_file = props.getString("density");
        heat_file = props.getString("heat");
        flame_file = props.getString("flame");
        read_data();

        VolumeHeteroFlameData flameData;
        flameData.toWorld = toWorld;
        flameData.toLocal = toLocal;
        flameData.density = density;
        flameData.density_index = density_index;
        flameData.flame = flame;
        flameData.flame_index = flame_index;
        flameData.heat = heat;
        flameData.heat_index = heat_index;
        VolumeEmitterHetero* tmp_emitter = new VolumeEmitterHetero(props);
        tmp_emitter->transfer_data(flameData);
        m_emitter = (Emitter *)tmp_emitter;
    }

    ~VolumeMediaHetero() {
        delete m_bsdf;
        release_data(density, density_index);
        release_data(flame, flame_index);
        release_data(heat, heat_index);
    }

    void read_data();
    void read_data(std::string filename, float*** &data, int* &data_index);
    void release_data(float*** &data, int* &data_index);

    std::string getName() const{
        return m_name;
    }

    float sample_dis(const Ray3f &ray, float t, Sampler* sampler) const {
        float total_dis = 0;
        float random, threshold, flame_value, density_value;
        Point3f pos;
        if (theta_t > Epsilon) {
            while (true) {
                float dis = sampler->next1D();
                dis = -(log(dis) / theta_t);
                total_dis += dis;
                if (total_dis >= t)
                    return std::numeric_limits<float>::infinity();
                pos = ray.o + total_dis * ray.d;
                pos = toLocal * pos;
                if (pos.norm() > 10)
                    return std::numeric_limits<float>::infinity();
                density_value = interpolation_3d(density, density_index, pos);
                if (density_value < -1.0f)
                    continue;
                flame_value = interpolation_3d(flame, flame_index, pos);
                if (flame_value < -1.0f)
                    continue;
                threshold = flame_value + density_value;
                random = sampler->next1D();
                if (random < threshold)
                    return total_dis;
            }
        }
        else
            return std::numeric_limits<float>::infinity();
    }

    BSDF* getBSDF() const {
        return m_bsdf;
    }

    Emitter* getEmitter() const {
        return m_emitter;
    }

    std::string toString() const {
        return tfm::format(
            "VolumeMedia Heterogeneous[\n"
            " theta_t = %f\n"
            "]",
            theta_t
        );
    }

private:
    std::string m_name;
    BSDF* m_bsdf = nullptr;
    float theta_t;
    Transform toWorld, toLocal;
    float*** density = nullptr;
    int* density_index = nullptr;
    float*** flame = nullptr;
    int* flame_index = nullptr;
    float*** heat = nullptr;
    int* heat_index = nullptr;
    std::string density_file;
    std::string flame_file;
    std::string heat_file;
    Emitter* m_emitter = nullptr;
};

void VolumeMediaHetero::read_data(std::string filename, float*** &data, int* &data_index) {
    std::ifstream fin(filename);
    data_index = new int[3];
    for (int i=0; i<3; i++) {
        fin >> data_index[i];
    }
    float value;
    data = new float**[data_index[0]];
    for (int i=0; i<data_index[0]; i++) {
        data[i] = new float*[data_index[1]];
        for (int j=0; j<data_index[1]; j++) {
            data[i][j] = new float[data_index[2]];
            for (int k=0; k<data_index[2]; k++) {
                fin >> value;
                data[i][j][k] = value;
            }
        }
    }
    fin.close();
}

void VolumeMediaHetero::read_data() {
    read_data(density_file, density, density_index);
    read_data(flame_file, flame, flame_index);
    read_data(heat_file, heat, heat_index);
}

void VolumeMediaHetero::release_data(float*** &data, int* &data_index) {
    if (data) {
        for (int i=0; i<data_index[0]; i++) {
            for (int j=0; j<data_index[1]; j++)
                delete data[i][j];
            delete data[i];
        }
        delete data;
        delete data_index;
    }
}


NORI_REGISTER_CLASS(VolumeMediaHetero, "volumemediahetero");
NORI_NAMESPACE_END
