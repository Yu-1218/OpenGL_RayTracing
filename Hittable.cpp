#include "Hittable.h"

// Sphere
bool Sphere::Hit(const Ray &ray, HitRecord *hit_record) const
{
    // TODO: Add your code here.
    bool ret = false;

    // Calculate the parameters for the equation
    float A = 1;
    float B = 2 * glm::dot((ray.o - this->o_), ray.d);
    float C = glm::dot(ray.o - this->o_, ray.o - this->o_) - pow(this->r_, 2);

    // std::cout << "o point " << o_[0] << " " << o_[1] << " " << o_[2] << std::endl;
    // std::cout << "radius " << r_ << std::endl;

    float delta = pow(B, 2) - 4 * A * C;
    if (delta >= 0)
    {
        ret = true;
    }
    else
    {
        ret = false;
    }

    if (float(((-1 * B - pow(delta, 0.5)) / (2 * A))) < 0)
    {
        ret = false;
    }

    if (ret)
    {
        hit_record->position = ray.o + float((-B - pow(delta, 0.5)) / (2 * A)) * ray.d;
        hit_record->normal = glm::normalize(hit_record->position - this->o_);
        hit_record->distance = glm::distance(ray.o, hit_record->position);
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(-2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal + ray.d);
        hit_record->material = material_;
    }
    return ret;
}

// Quadric
bool Quadric::Hit(const Ray &ray, HitRecord *hit_record) const
{
    // TODO: Add your code here.
    bool ret = false;

    // Calculate the parameters for the equation
    glm::vec4 O = glm::vec4(ray.o[0], ray.o[1], ray.o[2], 1);
    glm::vec4 D = glm::vec4(ray.d[0], ray.d[1], ray.d[2], 0);
    float A = glm::dot(D, this->A_ * D);
    float B = 2 * glm::dot(O, this->A_ * D);
    float C = glm::dot(O, this->A_ * O);

    float delta = pow(B, 2) - 4 * A * C;
    if (delta >= 0)
    {
        ret = true;
    }
    else
    {
        ret = false;
    }

    if (float(((-1 * B - pow(delta, 0.5)) / (2 * A))) < 0)
    {
        ret = false;
    }

    if (ret)
    {
        hit_record->position = ray.o + float(((-1 * B - pow(delta, 0.5)) / (2 * A))) * ray.d;
        hit_record->normal = (this->A_ + glm::transpose(this->A_)) * glm::vec4(hit_record->position[0], hit_record->position[1], hit_record->position[2], 1);
        hit_record->distance = glm::distance(ray.o, hit_record->position);
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(-2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal + ray.d);
        hit_record->material = material_;
    }
    return ret;
}

// Triangle
bool Triangle::Hit(const Ray &ray, HitRecord *hit_record) const
{
    // TODO: Add your code here.
    bool ret = false;

    // Calculate the parameters for the equation
    // N is the normal vector of the plane containing the triangle
    Vec N = glm::cross(this->b_ - this->a_, this->c_ - this->a_);
    float t = (glm::dot(this->c_ - ray.o, N)) / (glm::dot(ray.d, N));

    if (t < 0)
        return false;

    Vec o = ray.o + t * ray.d;
    Vec res_1 = glm::cross(this->a_ - o, this->b_ - o);
    Vec res_2 = glm::cross(this->b_ - o, this->c_ - o);
    Vec res_3 = glm::cross(this->c_ - o, this->a_ - o);
    if (glm::dot(res_1, res_2) > 0 && glm::dot(res_2, res_3) > 0 && glm::dot(res_1, res_3) > 0)
    {
        ret = true;
    }

    if (ret)
    {
        hit_record->position = ray.o + t * ray.d;
        if (phong_interpolation_)
        {
            // barycentric coordinates (u, v, w) for triangle (a, b, c)
            Vec v0 = this->b_ - this->a_;
            Vec v1 = this->c_ - this->a_;
            Vec v2 = o - this->a_;
            float d00 = glm::dot(v0, v0);
            float d01 = glm::dot(v0, v1);
            float d11 = glm::dot(v1, v1);
            float d20 = glm::dot(v2, v0);
            float d21 = glm::dot(v2, v1);
            float denom = d00 * d11 - d01 * d01;
            float v = (d11 * d20 - d01 * d21) / denom;
            float w = (d00 * d21 - d01 * d20) / denom;
            float u = 1.0f - v - w;
            hit_record->normal = u * this->n_a_ + v * this->n_b_ + w * this->n_c_;
        }
        else
        {
            hit_record->normal = glm::cross(this->b_ - this->a_, this->c_ - this->a_);
        }
        // no need to set material in this function
        hit_record->distance = glm::distance(ray.o, hit_record->position);
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(-2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal + ray.d);
    }
    return ret;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray &ray, HitRecord *hit_record) const
{
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret)
    {
        hit_record->material = material_;
    }
    return ret;
}

// Mesh
Mesh::Mesh(const std::string &file_path,
           const Material &material,
           bool phong_interpolation) : ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation)
{
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++)
    {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto &face : f_ind_)
    {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto &vertex_normal : vertex_normals_)
    {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto &face : f_ind_)
    {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min(1e5f, 1e5f, 1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto &vertex : vertices_)
    {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++)
    {
        InsertFace(root_, i);
    }
}

bool Mesh::Hit(const Ray &ray, HitRecord *hit_record) const
{
    const bool brute_force = false;
    if (brute_force)
    {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_)
        {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record))
            {
                if (curr_hit_record.distance < min_dist)
                {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f)
        {
            hit_record->material = material_;
            return true;
        }
        return false;
    }
    else
    {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret)
        {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t> &face, const Point &bbox_min, const Point &bbox_max) const
{
    for (size_t idx : face)
    {
        const auto &pt = vertices_[idx];
        for (int i = 0; i < 3; i++)
        {
            if (pt[i] < bbox_min[i] + 1e-6f)
                return false;
            if (pt[i] > bbox_max[i] - 1e-6f)
                return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray &ray, const Point &bbox_min, const Point &bbox_max) const
{
    float t_min = -1e5f;
    float t_max = 1e5f;

    for (int i = 0; i < 3; i++)
    {
        if (glm::abs(ray.d[i]) < 1e-6f)
        {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f)
            {
                t_min = 1e5f;
                t_max = -1e5f;
            }
        }
        else
        {
            if (ray.d[i] > 0.f)
            {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else
            {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode *u, size_t face_idx)
{
    const Point &bbox_min = u->bbox_min;
    const Point &bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++)
    {
        for (size_t b = 0; b < 2; b++)
        {
            for (size_t c = 0; c < 2; c++)
            {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max))
                {
                    if (u->childs[child_idx] == nullptr)
                    {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode *child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs)
    {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode *u, const Ray &ray, HitRecord *hit_record) const
{
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max))
    {
        return false;
    }
    float distance = 1e5f;
    for (const auto &face_idx : u->face_index)
    {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < distance)
            {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto &child : u->childs)
    {
        if (child == nullptr)
        {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < distance)
            {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}

// Hittable list
void HittableList::PushHittable(const Hittable &hittable)
{
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray &ray, HitRecord *hit_record) const
{
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_)
    {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < min_dist)
            {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}