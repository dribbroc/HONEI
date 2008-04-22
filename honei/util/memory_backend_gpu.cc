/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/memory_manager.hh>
#include <honei/util/memory_backend_gpu.hh>
#include <honei/util/tags.hh>

#include <GL/glew.h>
#include <GL/glxew.h>
#include <X11/X.h>
#include <X11/Xlib.h>

#include <cmath>

using namespace honei;

template class InstantiationPolicy<GPUBackend, Singleton>;

GPUBackend::Chunk::Chunk(unsigned long our_size) :
    size(our_size),
    width(GLuint(sqrt(double(our_size)))),
    fbo(0)
{
    if (size & (size - 1) != 0)
        throw MemoryChunkSizeInvalid(size, "Not a power of 2");

    glGenTextures(1, &id);
    glBindTexture(texture_target, id);
    glTexParameteri(texture_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(texture_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(texture_target, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(texture_target, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glTexImage2D(texture_target, 0, internal_format, width, width, 0,
            format, GL_FLOAT, 0);
}

GPUBackend::Chunk::Chunk(const GPUBackend::Chunk & other) :
    size(other.size),
    width(other.width),
    id(other.id)
{
}

GPUBackend::Chunk::~Chunk()
{
    glDeleteTextures(1, &id);
}

class GPUBackend::Data
{
    public:
        /// Our framebuffer object id. \todo Move to FramebufferObjectManager/Chunk
        GLuint fbo_id;

        /// Our set of known memory ids.
        std::set<MemoryId> known_ids;

        /// Our map of memory ids to GPU memory chunks.
        std::map<MemoryId, GPUBackend::Chunk *> chunk_map;

        /// Our set of known framebuffer objects.
        std::set<GLuint> fbos;

        /// Our map of framebuffer objects to bind-count/size.
        typedef std::map<GLuint, std::pair<unsigned short, unsigned long> > BindMap;
        BindMap bind_map;

        /// Our display.
        Display * display;

        /// Our window.
        Window window;

        /// Our OpenGL context.
        GLXContext gl_context;
};

GPUBackend::GPUBackend() :
    _data(new GPUBackend::Data)
{
    GLint attributes[2] = { GLX_RGBA, None };
    XSetWindowAttributes window_attributes;
    Window root_window;
    XVisualInfo* visual_info;
    int error(0);

    _data->display = XOpenDisplay(0);
    if(! _data->display)
        throw InternalError("Cannot connect to X server");

    root_window = DefaultRootWindow(_data->display);
    visual_info = glXChooseVisual(_data->display, 0, attributes);
    if (! visual_info)
        throw InternalError("No appropriate visual found");

    window_attributes.colormap = XCreateColormap(_data->display, root_window, visual_info->visual, AllocNone);
    _data->window = XCreateWindow(_data->display, root_window, 0, 0, 600, 600, 0, visual_info->depth, InputOutput,
            visual_info->visual, CWColormap, &window_attributes);
    _data->gl_context = glXCreateContext(_data->display, visual_info, 0, GL_TRUE);
    glXMakeCurrent(_data->display, _data->window, _data->gl_context);

    if (GLEW_OK != (error = glewInit()))
        throw InternalError("glewInit() returned '"
                + stringify(reinterpret_cast<const char *>(glewGetErrorString(error))) + "'");

    if (! glewIsSupported("GL_EXT_framebuffer_object"))
        throw InternalError("GL_EXT_framebuffer_object is missing");

    if (! glewIsSupported("GL_NV_float_buffer"))
        throw InternalError("GL_NV_float_buffer is missing");

    /// \todo Check for OpenGL 2.0 support

    glGenFramebuffersEXT(1, &_data->fbo_id);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _data->fbo_id);
}

GPUBackend::~GPUBackend()
{
    for (std::map<MemoryId, GPUBackend::Chunk *>::iterator c(_data->chunk_map.begin()),
            c_end(_data->chunk_map.end()) ; c != c_end ; ++c)
    {
        free(c->second);
    }

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glDeleteFramebuffersEXT(1, &_data->fbo_id);

    glXMakeCurrent(_data->display, None, 0);
    glXDestroyContext(_data->display, _data->gl_context);
    XDestroyWindow(_data->display, _data->window);
    XCloseDisplay(_data->display);
}

void
GPUBackend::_bind_framebuffer_object(Chunk * chunk)
{
    if (0 != chunk->fbo)
    {
        // We already have been bound to an fbo.
        /// \todo Throw an exception?
        return;
    }

    for (Data::BindMap::iterator i(_data->bind_map.begin()), i_end(_data->bind_map.end()) ;
            i != i_end ; ++i)
    {
        if (4 == i->second.first)
            continue;

        if (chunk->size != i->second.second)
            continue;

        ++(i->second.first);
        chunk->fbo = i->first;

        return;
    }

    /// \todo Allocate a new fbo.
}

void
GPUBackend::_release_framebuffer_object(Chunk * chunk)
{
    /// \todo implement
}

const GPUBackend::Chunk &
GPUBackend::get_chunk(const MemoryId id, const DeviceId)
{
    if (_data->known_ids.end() == _data->known_ids.find(id))
        throw MemoryIdNotKnown(id);

    std::map<MemoryId, GPUBackend::Chunk *>::const_iterator i(_data->chunk_map.find(id));
    ASSERT(_data->chunk_map.end() != i, "No information found for id '" + stringify(id) + "'!");

    return *(i->second);
}

MemoryBackend *
GPUBackend::backend_instance()
{
    return instance();
}

void
GPUBackend::upload(const MemoryId, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    /// \todo
}

void
GPUBackend::download(const MemoryId, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    /// \todo
}

void
GPUBackend::swap(const MemoryId left, const MemoryId right)
{
    /// \todo
}

void
GPUBackend::free(const MemoryId id, const DeviceId)
{
    std::map<MemoryId, GPUBackend::Chunk *>::iterator c(_data->chunk_map.find(id));
    if (_data->chunk_map.end() == c)
        free(c->second);
}

const GPUBackend::Chunk *
GPUBackend::alloc(unsigned long size)
{
    GPUBackend::Chunk * result(new GPUBackend::Chunk(size));

    // Do not record chunks in alloc().
    return result;
}

inline void
GPUBackend::free(const GPUBackend::Chunk * chunk)
{
    delete chunk;
}

static MemoryBackendRegistrator gpu_backend_registrator(tags::tv_gpu, &GPUBackend::backend_instance);
