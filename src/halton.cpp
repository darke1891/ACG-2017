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

#include <nori/sampler.h>
#include <nori/block.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

class Halton : public Sampler {
public:
    Halton(const PropertyList &propList) {
        m_sampleCount = (size_t) propList.getInteger("sampleCount", 1);
        seed_num = 4;
        scale = 1.0f - 2 * Epsilon;
        array_init();
    }

    virtual ~Halton() { 
        delete seeds;
        delete ran_3;
        delete ran_5;
    }

    void array_init() {
        seeds = new uint32_t[seed_num];
        for (int i=0; i<seed_num; i++)
            seeds[i] = m_random.nextUInt();
        seed_states = new uint32_t[seed_num];
        for (int i=0; i<seed_num; i++)
            seed_states[i] = 0;

        ran_3 = new uint32_t[3];
        ran_5 = new uint32_t[5];
        ran_7 = new uint32_t[7];
        ran_3[0] = 1;
        ran_3[1] = 2;
        ran_3[2] = 0;
        
        ran_5[0] = 1;
        ran_5[1] = 3;
        ran_5[2] = 4;
        ran_5[3] = 0;
        ran_5[4] = 2;

        ran_7[0] = 5;
        ran_7[1] = 6;
        ran_7[2] = 2;
        ran_7[3] = 4;
        ran_7[4] = 1;
        ran_7[5] = 0;
        ran_7[6] = 3;
    }

    std::unique_ptr<Sampler> clone() const {
        std::unique_ptr<Halton> cloned(new Halton());
        cloned->m_sampleCount = m_sampleCount;
        cloned->m_random = m_random;
        cloned->seed_num = seed_num;
        cloned->scale = scale;
        cloned->array_init();
        return std::move(cloned);
    }

    void prepare(const ImageBlock &block) {
        m_random.seed(
            block.getOffset().x(),
            block.getOffset().y()
        );
        for (int i=0; i<seed_num; i++)
            seeds[i] = m_random.nextUInt();
    }

    void generate() {
        for (int i=0; i<seed_num; i++)
            seeds[i] = m_random.nextUInt();
    }

    void advance()  { /* No-op for this sampler */ }

    float next1D() {
        return m_random.nextFloat();
    }
    
    Point2f next2D() {
        return Point2f(
            m_random.nextFloat(),
            m_random.nextFloat()
        );
    }

    float convert2(uint32_t n) {
        n = (n << 16) | (n >> 16);
        n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
        n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
        n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
        n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
        n = n ^ 0xffffffff;
        float res = n * 0x1p-32;
        return res;
    }

    float convert3(uint32_t n) {
        uint64_t n2 = n;
        uint64_t tmpn, tmpn2;
        uint64_t sbase = 1;
        uint64_t result = 0;
        uint64_t magic_number = 2863311531;
        for (int i=0; i<20; i++) {
            tmpn = (n2 * magic_number) >> 33;
            while (tmpn * 3 > n2)
                tmpn--;
            tmpn2 = n2 - tmpn * 3;
            tmpn2 = ran_3[tmpn2];
            result = result * 3 + 2 - tmpn2;
            sbase *= 3;
            n2 = tmpn;
        }
        float res = result;
        res /= sbase;
        return res;
    }

    float convert5(uint32_t n) {
        uint64_t n2 = n;
        uint64_t tmpn, tmpn2;
        uint64_t sbase = 1;
        uint64_t result = 0;
        uint64_t magic_number = 3435973837;
        for (int i=0; i<14; i++) {
            tmpn = (n2 * magic_number) >> 34;
            while (tmpn * 5 > n2)
                tmpn--;
            tmpn2 = n2 - tmpn * 5;
            tmpn2 = ran_5[tmpn2];
            result = result * 5 + tmpn2;
            sbase *= 5;
            n2 = tmpn;
        }
        float res = result;
        res /= sbase;
        return res;
    }

    float convert7(uint32_t n) {
        uint64_t n2 = n;
        uint64_t tmpn, tmpn2;
        uint64_t sbase = 1;
        uint64_t result = 0;
        uint64_t magic_number = 2454267027;
        for (int i=0; i<12; i++) {
            tmpn = (n2 * magic_number) >> 34;
            while (tmpn * 7 > n2)
                tmpn--;
            tmpn2 = n2 - tmpn * 7;
            tmpn2 = ran_7[tmpn2];
            result = result * 7 + tmpn2;
            sbase *= 7;
            n2 = tmpn;
        }
        float res = result;
        res /= sbase;
        return res;
    }

    Point2f next2D(int index) {
        if (index < seed_num) {
            uint32_t this_seed = seeds[index];
            seeds[index] += 1;
            switch(seed_states[index]) {
                case 0:
                    seed_states[index] = 1;
                    return Point2f(
                        convert2(this_seed) * scale + Epsilon,
                        convert3(this_seed) * scale + Epsilon
                    );
                case 1:
                default:
                    seed_states[index] = 0;
                    seeds[index] += 1;
                    return Point2f(
                        convert5(this_seed) * scale + Epsilon,
                        convert7(this_seed) * scale + Epsilon
                        );
            }
        }
        else
            return next2D();
    }

    std::string toString() const {
        return tfm::format("Halton[sampleCount=%i]", m_sampleCount);
    }
protected:
    Halton() { }

private:
    pcg32 m_random;
    uint32_t* seeds;
    uint32_t* seed_states;
    int seed_num;
    float scale;
    uint32_t *ran_3, *ran_5, *ran_7;
};

NORI_REGISTER_CLASS(Halton, "halton");
NORI_NAMESPACE_END
